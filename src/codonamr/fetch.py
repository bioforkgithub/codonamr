"""Sequence retrieval from NCBI, with an on-disk cache and a polite rate limit.

The product promise of this package is that the user supplies identifiers and
nothing else, so everything the analysis needs has to be discoverable from an
accession: the coding sequence itself, the source organism, and that organism's
genome, which supplies the host codon table the gene is compared against.

Network conduct follows the NCBI E-utilities usage policy: an email address is
mandatory, an API key is optional, and requests are limited to 3 per second
without a key and 10 per second with one. Every response is cached on disk under
the requesting accession, so a rerun of the same analysis makes no network calls
at all. Nothing here writes outside the cache directory.

References
----------
Sayers EW et al. (2022) Database resources of the NCBI. Nucleic Acids Res 50:D20-D26.
NCBI E-utilities usage guidelines, https://www.ncbi.nlm.nih.gov/books/NBK25497/
"""
import gzip
import io
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass, field
from pathlib import Path

from Bio import SeqIO

__all__ = ["CdsRecord", "fetch_cds", "read_fasta", "resolve_host_genome",
           "default_cache_dir", "guess_db", "AssemblyHit"]

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq"
USER_AGENT = "codonamr/0.1.0 (https://github.com/; polite E-utilities client)"

#: how many accessions go into one efetch call
BATCH_SIZE = 200

#: protein accession prefixes that have to be resolved back to a nucleotide CDS
_PROTEIN_PREFIXES = ("WP_", "NP_", "YP_", "XP_", "AP_", "ZP_")
_PROTEIN_RE = re.compile(r"^[A-Z]{3}[0-9]{5,7}(\.[0-9]+)?$")
_NUCL_RE = re.compile(r"^[A-Z]{1,2}[0-9]{5,8}(\.[0-9]+)?$")


# ---------------------------------------------------------------- records
@dataclass
class CdsRecord:
    """One coding sequence and the provenance needed to interpret it.

    ``sequence`` is raw, as deposited. It has not been through QC yet: that is
    :func:`codonamr.qc.check_cds`, kept separate so rejection reasons are
    reported rather than applied silently during download.
    """
    id: str
    sequence: str
    description: str = ""
    organism: str = None
    taxid: str = None
    product: str = None
    gene: str = None
    source: str = "ncbi"
    extra: dict = field(default_factory=dict)

    def as_dict(self):
        return asdict(self)


# ------------------------------------------------------------------ cache
def default_cache_dir():
    """Cache location, honouring ``CODONAMR_CACHE`` then ``XDG_CACHE_HOME``.

    Nothing is hardcoded. If neither variable is set this falls back to the
    per-user cache directory the platform already uses for such things.
    """
    env = os.environ.get("CODONAMR_CACHE")
    if env:
        return Path(env).expanduser()
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg).expanduser() if xdg else Path.home() / ".cache"
    return base / "codonamr"


def _cache_file(cache_dir, namespace, key, suffix=""):
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", str(key))[:150]
    d = Path(cache_dir) / namespace
    d.mkdir(parents=True, exist_ok=True)
    return d / (safe + suffix)


def _read_cache(path, binary=False):
    try:
        return path.read_bytes() if binary else path.read_text()
    except OSError:
        return None


def _write_cache(path, data):
    tmp = path.with_suffix(path.suffix + ".part")
    if isinstance(data, bytes):
        tmp.write_bytes(data)
    else:
        tmp.write_text(data)
    tmp.replace(path)


# ------------------------------------------------------------------- http
class _Throttle:
    """Minimum spacing between requests, shared by every call in a process."""

    def __init__(self, per_second):
        self.interval = 1.0 / float(per_second)
        self.last = 0.0

    def wait(self):
        gap = self.interval - (time.monotonic() - self.last)
        if gap > 0:
            time.sleep(gap)
        self.last = time.monotonic()


_THROTTLES = {}


def _throttle_for(api_key):
    """3 requests/second without a key, 10 with one, per the NCBI policy."""
    rate = 10.0 if api_key else 3.0
    if rate not in _THROTTLES:
        _THROTTLES[rate] = _Throttle(rate)
    return _THROTTLES[rate]


def _http_get(url, params=None, api_key=None, retries=5, timeout=120,
              binary=False, verbose=True):
    """GET with backoff on 429 and 5xx, returning text or bytes.

    Backoff is exponential starting at 1 second. HTTP 429 is what NCBI returns
    when the rate limit is exceeded, so it is retried rather than raised, but
    the throttle above should keep it from happening in the first place.
    """
    if params:
        url = url + ("&" if "?" in url else "?") + urllib.parse.urlencode(params)
    _throttle_for(api_key).wait()
    last = None
    for attempt in range(retries):
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                raw = resp.read()
            return raw if binary else raw.decode("utf-8", "replace")
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code not in (429, 500, 502, 503, 504):
                raise
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
        delay = 2.0 ** attempt
        if verbose:
            print("codonamr: retrying in %.0fs after %s" % (delay, last),
                  file=sys.stderr)
        time.sleep(delay)
    raise RuntimeError("giving up on %s after %d attempts: %s"
                       % (url.split("?")[0], retries, last))


def _require_email(email):
    if not email:
        raise ValueError(
            "NCBI requires an email address with every E-utilities request. "
            "Pass email=... , set the CODONAMR_EMAIL environment variable, or "
            "supply local sequences with read_fasta() instead."
        )
    return email


def _eutils(tool, email, api_key=None, **params):
    params.update(tool="codonamr", email=_require_email(email))
    if api_key:
        params["api_key"] = api_key
    return _http_get(EUTILS + tool + ".fcgi", params, api_key=api_key)


# ------------------------------------------------------------ id handling
def guess_db(accession):
    """Guess ``"protein"`` or ``"nuccore"`` from the shape of an accession.

    RefSeq prefixes are unambiguous. For INSDC accessions the rule used here is
    that three letters followed by digits is a protein and one or two letters
    followed by digits is a nucleotide. That is the usual convention but it is a
    convention, not a guarantee, so pass ``db=`` explicitly if it matters.
    """
    a = str(accession).strip().upper()
    if a.startswith(_PROTEIN_PREFIXES):
        return "protein"
    if a[:3] in ("NC_", "NZ_", "NM_", "NR_", "NG_", "NT_", "NW_", "AC_"):
        return "nuccore"
    if _PROTEIN_RE.match(a):
        return "protein"
    return "nuccore"


def _parse_ids(ids):
    out = []
    for i in ids:
        i = str(i).strip()
        if i and not i.startswith("#"):
            out.append(i)
    return out


# ------------------------------------------------------- GenBank handling
def _split_genbank(text):
    """Split a multi-record GenBank or GenPept stream into (accession, text)."""
    out, buf = [], []
    for line in text.splitlines(keepends=True):
        buf.append(line)
        if line.rstrip() == "//":
            chunk = "".join(buf)
            buf = []
            acc = None
            for l in chunk.splitlines():
                if l.startswith("VERSION"):
                    acc = l.split()[1]
                    break
                if l.startswith("ACCESSION"):
                    acc = l.split()[1]
            if acc:
                out.append((acc, chunk))
    return out


def _source_taxid(rec):
    for f in rec.features:
        if f.type == "source":
            for x in f.qualifiers.get("db_xref", []):
                if x.lower().startswith("taxon:"):
                    return x.split(":", 1)[1]
    return None


def _first(q, key):
    v = q.get(key)
    return v[0] if v else None


def _records_from_genbank(text, requested_id):
    """Turn cached GenBank text into :class:`CdsRecord` objects.

    A record with one annotated CDS yields one sequence. A record with several,
    which happens when an accession points at a plasmid or a whole genome rather
    than at a gene, yields one record per CDS with a suffixed id. A record with
    no CDS annotation yields its whole sequence, since that is usually a
    submitted gene-only entry, and QC will catch it if it is not in frame.
    """
    out = []
    for rec in SeqIO.parse(io.StringIO(text), "genbank"):
        organism = rec.annotations.get("organism")
        taxid = _source_taxid(rec)
        try:
            full = str(rec.seq)
        except Exception:
            full = ""
        cds = [f for f in rec.features if f.type == "CDS"]
        if cds and full:
            for i, f in enumerate(cds, 1):
                if "pseudo" in f.qualifiers or "pseudogene" in f.qualifiers:
                    continue
                try:
                    sub = str(f.extract(rec.seq))
                except Exception:
                    continue
                pid = _first(f.qualifiers, "protein_id")
                sid = rec.id if len(cds) == 1 else "%s_cds_%s" % (rec.id, pid or i)
                out.append(CdsRecord(
                    id=sid, sequence=sub, description=rec.description,
                    organism=organism, taxid=taxid,
                    product=_first(f.qualifiers, "product"),
                    gene=_first(f.qualifiers, "gene"),
                    extra={"requested": requested_id, "parent": rec.id}))
        elif full:
            out.append(CdsRecord(
                id=rec.id, sequence=full, description=rec.description,
                organism=organism, taxid=taxid,
                extra={"requested": requested_id, "parent": rec.id,
                       "note": "no CDS feature, whole record used"}))
    return out


# --------------------------------------------------- protein to nucleotide
_CODED_BY = re.compile(
    r"(?P<comp>complement\()?(?P<acc>[A-Za-z0-9_.]+):"
    r"(?P<start><?\d+)\.\.(?P<stop>>?\d+)")


def _coded_by_regions(text):
    """Parse ``/coded_by=`` from a GenPept record into (acc, start, stop, comp).

    ``join(...)`` is handled by taking every listed region in order, which is
    what is wanted for a programmed frameshift or a spliced eukaryotic CDS.
    Fuzzy ends (``<`` or ``>``) are kept but flagged, because a partial CDS is a
    QC failure rather than something to silently trim.
    """
    m = re.search(r'/coded_by="([^"]+)"', text.replace("\n", " ").replace("  ", " "))
    if not m:
        m = re.search(r"coded_by=\"?([^\"\n]+)", text)
    if not m:
        return [], False
    spec = re.sub(r"\s+", "", m.group(1))
    regions, fuzzy = [], ("<" in spec or ">" in spec)
    for g in _CODED_BY.finditer(spec):
        regions.append((g.group("acc"), int(g.group("start").lstrip("<")),
                        int(g.group("stop").lstrip(">")), bool(g.group("comp"))))
    return regions, fuzzy


def _fetch_region(acc, start, stop, complement, email, api_key, cache_dir):
    key = "%s_%d_%d_%d" % (acc, start, stop, int(complement))
    path = _cache_file(cache_dir, "region", key, ".fasta")
    text = _read_cache(path)
    if text is None:
        text = _eutils("efetch", email, api_key, db="nuccore", id=acc,
                       rettype="fasta", retmode="text",
                       seq_start=start, seq_stop=stop,
                       strand=2 if complement else 1)
        _write_cache(path, text)
    recs = list(SeqIO.parse(io.StringIO(text), "fasta"))
    return str(recs[0].seq) if recs else ""


def _resolve_protein(acc, text, email, api_key, cache_dir):
    regions, fuzzy = _coded_by_regions(text)
    if not regions:
        return None
    parts = [_fetch_region(a, s, e, c, email, api_key, cache_dir)
             for a, s, e, c in regions]
    seq = "".join(p for p in parts if p)
    if not seq:
        return None
    rec = next(SeqIO.parse(io.StringIO(text), "genbank"), None)
    organism = rec.annotations.get("organism") if rec else None
    taxid = _source_taxid(rec) if rec else None
    product = gene = None
    if rec:
        for f in rec.features:
            if f.type == "CDS":
                product = _first(f.qualifiers, "product") or product
                gene = _first(f.qualifiers, "gene") or gene
    return CdsRecord(
        id=acc, sequence=seq, description=rec.description if rec else "",
        organism=organism, taxid=taxid, product=product, gene=gene,
        extra={"requested": acc, "resolved_from": "protein",
               "coded_by": regions[0][0], "partial": fuzzy})


# ------------------------------------------------------------- public API
def fetch_cds(ids, email, api_key=None, cache_dir=None, db="nuccore",
              batch_size=BATCH_SIZE, verbose=True):
    """Fetch coding sequences for a list of NCBI accessions.

    Parameters
    ----------
    ids : iterable of str
        Nucleotide accessions. Protein accessions are accepted too and are
        resolved to their coding nucleotide region through the ``/coded_by``
        qualifier of the protein record, which is exact, rather than through a
        link to the parent record, which is not.
    email : str
        Contact address sent with every request. Required by NCBI.
    api_key : str, optional
        NCBI API key. Raises the permitted rate from 3 to 10 requests/second.
    cache_dir : path-like, optional
        Where fetched records are stored. Defaults to :func:`default_cache_dir`.
    db : {"nuccore", "protein", "auto"}
        ``"auto"`` guesses per accession with :func:`guess_db`.
    batch_size : int
        Accessions per efetch call. 200 is the largest NCBI recommends for GET.

    Returns
    -------
    list of CdsRecord

    Notes
    -----
    Records are cached under the accession as requested, so rerunning an
    analysis fetches nothing. Delete the cache directory to force a refresh.
    One accession can yield several records if it carries several CDS features,
    so do not assume the output is the same length as the input.
    """
    cache_dir = Path(cache_dir) if cache_dir else default_cache_dir()
    wanted = _parse_ids(ids)
    if not wanted:
        return []
    email = _require_email(email)

    groups = {}
    for acc in wanted:
        target = guess_db(acc) if db == "auto" else db
        groups.setdefault(target, []).append(acc)

    out, missing = [], []
    cached = {}
    for target, accs in groups.items():
        todo = []
        for acc in accs:
            path = _cache_file(cache_dir, target, acc, ".gb")
            text = _read_cache(path)
            if text:
                cached[(target, acc)] = text
            else:
                todo.append(acc)
        for i in range(0, len(todo), batch_size):
            chunk = todo[i:i + batch_size]
            if verbose:
                print("codonamr: fetching %d %s records (%d/%d)"
                      % (len(chunk), target, i + len(chunk), len(todo)),
                      file=sys.stderr)
            text = _eutils("efetch", email, api_key, db=target,
                           id=",".join(chunk),
                           rettype="gp" if target == "protein" else "gb",
                           retmode="text")
            pieces = dict(_split_genbank(text))
            base = {k.split(".")[0]: v for k, v in pieces.items()}
            for acc in chunk:
                piece = pieces.get(acc) or base.get(acc.split(".")[0])
                if piece is None:
                    missing.append(acc)
                    continue
                _write_cache(_cache_file(cache_dir, target, acc, ".gb"), piece)
                cached[(target, acc)] = piece

    for target, accs in groups.items():
        for acc in accs:
            text = cached.get((target, acc))
            if text is None:
                if acc not in missing:
                    missing.append(acc)
                continue
            if target == "protein":
                rec = _resolve_protein(acc, text, email, api_key, cache_dir)
                if rec is None:
                    missing.append(acc)
                else:
                    out.append(rec)
            else:
                got = _records_from_genbank(text, acc)
                if not got:
                    missing.append(acc)
                out.extend(got)

    if missing and verbose:
        print("codonamr: %d accession(s) returned nothing: %s"
              % (len(missing), ", ".join(sorted(set(missing))[:10])),
              file=sys.stderr)
    return out


def read_fasta(path, organism=None, taxid=None):
    """Read local sequences, so the package is usable with no network at all.

    ``[organism=...]`` and ``[taxid=...]`` tags in the description are read if
    present, which is the format NCBI's own CDS downloads use. The ``organism``
    and ``taxid`` arguments override them for the whole file.
    """
    out = []
    for rec in SeqIO.parse(str(path), "fasta"):
        tags = dict(re.findall(r"\[(\w+)=([^\]]*)\]", rec.description))
        out.append(CdsRecord(
            id=rec.id, sequence=str(rec.seq).upper(),
            description=rec.description,
            organism=organism or tags.get("organism"),
            taxid=taxid or tags.get("taxid") or tags.get("db_xref", "").replace("taxon:", "") or None,
            product=tags.get("protein"), gene=tags.get("gene"),
            source="local"))
    return out


# ----------------------------------------------------------- host genomes
@dataclass
class AssemblyHit:
    """One row of an NCBI assembly summary file, plus the CDS it points at."""
    accession: str
    organism: str
    taxid: str
    species_taxid: str
    category: str
    level: str
    ftp_path: str
    division: str


_SUMMARY_COLUMNS = ("assembly_accession", "refseq_category", "taxid",
                    "species_taxid", "organism_name", "version_status",
                    "assembly_level", "ftp_path")
_SUMMARY_FALLBACK = {"assembly_accession": 0, "refseq_category": 4, "taxid": 5,
                     "species_taxid": 6, "organism_name": 7,
                     "version_status": 10, "assembly_level": 11, "ftp_path": 19}


def _summary_index(header_line):
    names = header_line.lstrip("#").rstrip("\n").split("\t")
    idx = {}
    for col in _SUMMARY_COLUMNS:
        idx[col] = names.index(col) if col in names else _SUMMARY_FALLBACK[col]
    return idx


def _assembly_table(division, cache_dir, max_age_days=30, verbose=True):
    """Download and cache the reference/representative rows for one division.

    The full RefSeq assembly summary for bacteria is hundreds of megabytes, so
    it is streamed and filtered line by line and only the reference and
    representative rows are kept on disk. That is a few thousand rows, which is
    what the host lookup actually needs.
    """
    path = _cache_file(cache_dir, "assembly", division, ".tsv")
    if path.exists():
        age = (time.time() - path.stat().st_mtime) / 86400.0
        if age < max_age_days:
            return [l.rstrip("\n").split("\t") for l in path.read_text().splitlines()]
    url = "%s/%s/assembly_summary.txt" % (FTP_BASE, division)
    if verbose:
        print("codonamr: downloading assembly summary for %s" % division,
              file=sys.stderr)
    _throttle_for(None).wait()
    rows, idx = [], None
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=600) as resp:
        for raw in io.TextIOWrapper(resp, encoding="utf-8", errors="replace"):
            if raw.startswith("#"):
                if "assembly_accession" in raw:
                    idx = _summary_index(raw)
                continue
            if idx is None:
                idx = dict(_SUMMARY_FALLBACK)
            f = raw.rstrip("\n").split("\t")
            if len(f) <= idx["ftp_path"]:
                continue
            if f[idx["refseq_category"]] not in ("reference genome",
                                                 "representative genome"):
                continue
            if f[idx["version_status"]] != "latest":
                continue
            rows.append([f[idx[c]] for c in _SUMMARY_COLUMNS])
    _write_cache(path, "\n".join("\t".join(r) for r in rows) + "\n")
    return rows


def _match_assembly(rows, division, taxid=None, organism=None):
    """Exact taxid, then species taxid, then organism name, then genus."""
    hits = []
    for r in rows:
        acc, cat, tid, sp, name, _status, level, ftp = r
        hit = AssemblyHit(acc, name, tid, sp, cat, level, ftp, division)
        if taxid and tid == str(taxid):
            hits.append((0, hit))
        elif taxid and sp == str(taxid):
            hits.append((1, hit))
        elif organism and name.lower() == organism.lower():
            hits.append((2, hit))
        elif organism and name.lower().startswith(organism.lower()):
            hits.append((3, hit))
        elif organism and organism.lower().split()[0] == name.lower().split()[0]:
            hits.append((4, hit))
    if not hits:
        return None
    hits.sort(key=lambda t: (t[0], t[1].category != "reference genome",
                             t[1].level != "Complete Genome"))
    return hits[0][1]


def _taxid_lineage(taxid, email, api_key, cache_dir):
    path = _cache_file(cache_dir, "taxonomy", taxid, ".json")
    text = _read_cache(path)
    if text is None:
        text = _eutils("esummary", email, api_key, db="taxonomy", id=str(taxid),
                       retmode="json")
        _write_cache(path, text)
    try:
        data = json.loads(text)
        rec = data["result"][str(taxid)]
        return rec.get("scientificname"), rec.get("division", "")
    except Exception:
        return None, ""


_DIVISION_MAP = {"bacteria": "bacteria", "archaea": "archaea",
                 "viruses": "viral", "phages": "viral", "fungi": "fungi",
                 "plants": "plant", "invertebrates": "invertebrate",
                 "rodents": "vertebrate_mammalian",
                 "mammals": "vertebrate_mammalian",
                 "primates": "vertebrate_mammalian",
                 "vertebrates": "vertebrate_other"}


def parse_cds_fasta(handle, organism=None, taxid=None, skip_pseudo=True):
    """Parse an NCBI ``_cds_from_genomic.fna`` stream into records.

    Headers carry bracketed tags such as ``[gene=rpsA]`` and
    ``[protein=30S ribosomal protein S1]``. The protein tag is what
    :func:`codonamr.host.reference_set` matches on, so it is kept verbatim.
    Pseudogenes are dropped by default: they are frameshifted or truncated by
    definition and would fail QC anyway, but dropping them here keeps the
    rejection report about the user's own genes.
    """
    out = []
    for rec in SeqIO.parse(handle, "fasta"):
        tags = dict(re.findall(r"\[(\w+)=([^\]]*)\]", rec.description))
        if skip_pseudo and tags.get("pseudo") == "true":
            continue
        out.append(CdsRecord(
            id=rec.id, sequence=str(rec.seq).upper(),
            description=rec.description,
            organism=tags.get("organism") or organism,
            taxid=taxid, product=tags.get("protein"), gene=tags.get("gene"),
            source="refseq_assembly",
            extra={"locus_tag": tags.get("locus_tag"),
                   "protein_id": tags.get("protein_id")}))
    return out


def resolve_host_genome(taxid_or_organism, email=None, api_key=None,
                        cache_dir=None, divisions=("bacteria", "archaea"),
                        max_cds=None, verbose=True):
    """Find the RefSeq reference or representative genome for an organism.

    This is what makes the host comparison automatic: given only the taxid or
    organism name attached to a fetched gene, it locates a genome for that
    organism and downloads its CDS set, which is then used to build the host
    codon table, the CAI reference set and the length-matched null.

    Parameters
    ----------
    taxid_or_organism : int or str
        A numeric NCBI taxid, or an organism name such as ``"Escherichia coli"``.
    divisions : tuple of str
        RefSeq genome divisions to search, in order. The default covers
        prokaryotes, which is where mobile resistance genes live. Add
        ``"viral"``, ``"fungi"`` or a vertebrate division for other hosts.
    max_cds : int, optional
        Truncate the returned CDS list. Only useful for tests.

    Returns
    -------
    dict
        ``{"assembly": AssemblyHit, "cds": [CdsRecord, ...]}``, or ``None`` if
        no reference or representative genome exists for the organism.

    Notes
    -----
    Matching falls back from exact taxid to species taxid to exact organism
    name to genus. A genus-level fallback means the host codon table comes from
    a relative, not the actual strain, which is usually acceptable for codon
    usage but is recorded in the returned ``AssemblyHit`` so the report can say
    so. Only reference and representative genomes are considered, which avoids
    picking one of the thousands of near-identical clinical isolates of a
    well-sequenced species, but also means an organism with no designated
    representative returns ``None`` rather than a silently arbitrary genome.
    """
    cache_dir = Path(cache_dir) if cache_dir else default_cache_dir()
    key = str(taxid_or_organism).strip()
    taxid = key if key.isdigit() else None
    organism = None if taxid else key
    if taxid and email:
        name, div = _taxid_lineage(taxid, email, api_key, cache_dir)
        organism = name or organism
        mapped = _DIVISION_MAP.get(str(div).lower())
        if mapped:
            divisions = (mapped,) + tuple(d for d in divisions if d != mapped)

    hit = None
    for division in divisions:
        try:
            rows = _assembly_table(division, cache_dir, verbose=verbose)
        except Exception as exc:
            if verbose:
                print("codonamr: could not read %s summary: %s" % (division, exc),
                      file=sys.stderr)
            continue
        hit = _match_assembly(rows, division, taxid, organism)
        if hit:
            break
    if hit is None:
        return None

    base = hit.ftp_path.replace("ftp://", "https://")
    url = "%s/%s_cds_from_genomic.fna.gz" % (base, base.rsplit("/", 1)[-1])
    path = _cache_file(cache_dir, "genome", hit.accession, "_cds.fna.gz")
    if not path.exists():
        if verbose:
            print("codonamr: downloading CDS set for %s (%s)"
                  % (hit.organism, hit.accession), file=sys.stderr)
        _write_cache(path, _http_get(url, binary=True, verbose=verbose))
    with gzip.open(path, "rt") as fh:
        cds = parse_cds_fasta(fh, organism=hit.organism, taxid=hit.taxid)
    if max_cds:
        cds = cds[:max_cds]
    return {"assembly": hit, "cds": cds}
