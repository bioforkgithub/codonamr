"""Host codon usage: reference sets, optimal codons, tAI and gene-to-host distance.

A codon-usage index is meaningless on its own. ENC says how biased a gene is,
but not whether it is biased in the direction its host prefers, and that is the
question a mobile resistance gene raises. Everything here builds the host side
of that comparison from the host's own genome.

References
----------
CAI reference sets  Sharp PM, Li WH (1987) Nucleic Acids Res 15:1281-1295.
Optimal codons      Ikemura T (1981) J Mol Biol 146:1-21.
tAI                 dos Reis M, Savva R, Wernisch L (2004) Nucleic Acids Res 32:5036-5044.
tRNAscan-SE         Chan PP, Lowe TM (2019) Methods Mol Biol 1962:1-14.
CUFS                Diament A, Pinter RY, Tuller T (2014) Nat Commun 5:5876.
Endres-Schindelin   Endres DM, Schindelin JE (2003) IEEE Trans Inf Theory 49:1858-1860.
"""
import math
import re
from collections import defaultdict

from .genetic_code import CODON2AA, CODONS, FAMILY, codon_list
from .metrics import build_w, codon_counts, rscu

__all__ = ["host_codon_table", "reference_set", "optimal_codons", "tai",
           "trna_copy_numbers_from_tRNAscan", "cufs", "rscu_distance",
           "TAI_S_VALUES"]


# ------------------------------------------------------------ host tables
def host_codon_table(cds_seqs):
    """Pooled codon counts, RSCU and relative adaptiveness for a genome.

    Parameters
    ----------
    cds_seqs : iterable of str
        QC-passed, in-frame coding sequences, normally a whole RefSeq CDS set.

    Returns
    -------
    dict
        ``counts``, ``rscu``, ``w``, ``n_genes`` and ``n_codons``.

    Notes
    -----
    The ``w`` returned here is genome-wide, not the CAI ``w``. CAI is defined
    against a set of highly expressed genes (see :func:`reference_set`), and a
    genome-wide ``w`` measures something different: adaptation to the bulk
    composition of the genome, which in a skewed genome is largely mutational
    rather than selective. It is provided as an explicit fallback for organisms
    where no reference set can be identified, and the report says when that
    fallback was used.
    """
    counts = {c: 0 for c in CODONS if CODON2AA[c] != "*"}
    n_genes = 0
    for s in cds_seqs:
        n_genes += 1
        for c, n in codon_counts(s).items():
            counts[c] += n
    return {"counts": counts, "rscu": rscu(counts), "w": build_w(cds_seqs),
            "n_genes": n_genes, "n_codons": sum(counts.values())}


# ------------------------------------------------------------ reference set
#: gene names that contain a ribosomal-protein name but encode an enzyme that
#: acts *on* a ribosomal protein. These are not highly expressed and including
#: them contaminates the reference set.
RP_EXCLUDED_GENES = frozenset("""
prmA prmB prmC rimI rimJ rimK rimL rimM rimO rimP roxA ycaO efp epmA epmB
rlmA rlmB rlmC rlmD rlmE rlmF rlmG rlmH rlmI rlmJ rlmKL rlmL rlmM rlmN
rsmA rsmB rsmC rsmD rsmE rsmF rsmG rsmH rsmI rsmJ rluA rluB rluC rluD rluE rluF
rbfA rsgA rsfS rimB ribosomal
""".split())

#: phrases that mark a product as a modifier, chaperone or assembly factor
#: rather than a structural ribosomal protein
RP_EXCLUDED_PHRASES = (
    "methyltransferase", "methylthiotransferase", "acetyltransferase",
    "transferase", "hydroxylase", "oxygenase", "ligase", "synthetase",
    "kinase", "phosphatase", "glycosylase", "deacetylase", "amidotransferase",
    "modification", "maturation", "biogenesis", "assembly", "chaperone",
    "pseudouridine", "rrna", "ribosome-binding", "ribosome binding",
    "recycling", "silencing", "hibernation", "gtpase", "atpase",
    "ribosomal protein serine acetyltransferase",
    "ribosomal-protein", "alanine n-acetyltransferase",
)

_RP_PRODUCT = re.compile(r"\b(30S|50S|40S|60S)?\s*ribosomal protein\b", re.I)
_RP_GENE = re.compile(r"^(rps|rpl|rpm|mrps|mrpl)[A-Z]\d*$", re.I)

#: optional additions to the reference set, following the classic E. coli set
ELONGATION_FACTOR_GENES = frozenset("tufA tufB tuf fusA tsf infA infB infC".split())


def _looks_like_ribosomal_protein(product, gene):
    g = (gene or "").strip()
    p = (product or "").strip()
    low = p.lower()
    if g and g.lower() in RP_EXCLUDED_GENES:
        return False
    if any(x in low for x in RP_EXCLUDED_PHRASES):
        return False
    if g and _RP_GENE.match(g):
        return True
    if _RP_PRODUCT.search(p):
        # "ribosomal protein L11 methyltransferase" is already excluded above;
        # what is left is a structural subunit name such as "50S ribosomal
        # protein L7/L12".
        return True
    return False


def reference_set(cds_seqs, products, genes=None, min_genes=20,
                  include_elongation_factors=False, extra_patterns=()):
    """Pick a highly expressed reference set for CAI, from product annotation.

    Parameters
    ----------
    cds_seqs : sequence of str
        Host coding sequences, already QC-passed.
    products : sequence of str or None
        Product annotation for each sequence, in the same order. These are the
        ``[protein=...]`` tags of an NCBI CDS download.
    genes : sequence of str or None, optional
        Gene symbols in the same order. Used when a product string is missing.
    min_genes : int
        Below this many hits the set is considered untrustworthy and an empty
        list is returned, so the caller can fall back to a genome-wide ``w``
        instead of computing a CAI from five genes.
    include_elongation_factors : bool
        Add translation elongation and initiation factors. Sharp and Li's
        original E. coli set combined ribosomal proteins, elongation factors
        and outer-membrane proteins. Ribosomal proteins alone are used by
        default because they are the part that transfers cleanly to an
        arbitrary prokaryote.
    extra_patterns : iterable of str
        Additional case-insensitive regular expressions matched against the
        product string.

    Returns
    -------
    list of str

    Notes
    -----
    The hard part is that ribosomal-protein *modification enzymes* carry the
    words "ribosomal protein" in their product names: prmA and prmB are L11 and
    L3 methyltransferases, rimK ligates glutamate onto S6, rimO is an S12
    methylthiotransferase, roxA hydroxylates L16 and ycaO is its paralogue.
    They are ordinary low-abundance enzymes. Letting them into the reference set
    pulls ``w`` towards average genome usage in exactly the families where the
    real ribosomal proteins are most biased, which flattens CAI for every gene
    scored afterwards. They are excluded by gene symbol and by product phrase,
    and the count of excluded near-misses is worth checking on an unfamiliar
    genome.
    """
    extra = [re.compile(p, re.I) for p in extra_patterns]
    genes = list(genes) if genes is not None else [None] * len(cds_seqs)
    out = []
    for seq, product, gene in zip(cds_seqs, products, genes):
        keep = _looks_like_ribosomal_protein(product, gene)
        if not keep and include_elongation_factors and gene:
            keep = gene in ELONGATION_FACTOR_GENES
        if not keep and extra and product:
            keep = any(r.search(product) for r in extra)
        if keep:
            out.append(seq)
    return out if len(out) >= min_genes else []


def optimal_codons(rscu_table, threshold=1.0):
    """Codons the host over-uses relative to even usage, that is RSCU > 1.

    This is the operational definition used for Fop when no expression data are
    available. Ikemura's original optimal codons were defined against measured
    tRNA abundances, and the RSCU rule is a proxy for that. In a genome with
    strong mutational GC skew the two definitions come apart, because a codon
    can be common without being translationally preferred.
    """
    return frozenset(c for c, v in rscu_table.items()
                     if v is not None and v > threshold
                     and len(FAMILY[CODON2AA[c]]) > 1)


# ---------------------------------------------------------------------- tAI
#: selective constraints of dos Reis, Savva & Wernisch (2004), Table 1. Keys are
#: (anticodon wobble base 34, codon base 3). A value of 0 means the pairing is
#: unpenalised; 0.9999 means it is possible but almost never used.
TAI_S_VALUES = {
    ("I", "U"): 0.0,     # A34 read as inosine, pairs with U3 as efficiently as A:U
    ("G", "C"): 0.0,
    ("U", "A"): 0.0,
    ("C", "G"): 0.0,
    ("G", "U"): 0.41,    # wobble G34:U3
    ("I", "C"): 0.28,    # wobble I34:C3
    ("I", "A"): 0.9999,  # wobble I34:A3, effectively no decoding
    ("U", "G"): 0.68,    # wobble U34:G3
    ("L", "A"): 0.89,    # lysidine-modified C34 of bacterial tRNA-Ile2 on AUA
}

_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A", "N": "N"}

#: key under which :func:`trna_copy_numbers_from_tRNAscan` stores the bacterial
#: lysidine tRNA-Ile, whose anticodon CAU is shared with tRNA-Met
ILE2_KEY = "CAT_Ile2"


def _revcomp(s):
    return "".join(_COMP.get(b, "N") for b in reversed(s.upper()))


def _tai_weights(tgcn, prokaryote=True):
    """Absolute adaptiveness W for every sense codon, before normalisation."""
    s = TAI_S_VALUES
    w = {}
    for c in CODONS:
        if CODON2AA[c] == "*":
            continue
        third = c[2]
        a_wc = _revcomp(c)                      # anticodon pairing Watson-Crick
        if third == "T":                        # codon U3: A34 (as I) and G34
            w[c] = ((1 - s[("I", "U")]) * tgcn.get(a_wc, 0)
                    + (1 - s[("G", "U")]) * tgcn.get(_revcomp(c[:2] + "C"), 0))
        elif third == "C":                      # codon C3: G34 and I34
            w[c] = ((1 - s[("G", "C")]) * tgcn.get(a_wc, 0)
                    + (1 - s[("I", "C")]) * tgcn.get(_revcomp(c[:2] + "T"), 0))
        elif third == "A":                      # codon A3: U34 and I34
            w[c] = ((1 - s[("U", "A")]) * tgcn.get(a_wc, 0)
                    + (1 - s[("I", "A")]) * tgcn.get(_revcomp(c[:2] + "T"), 0))
        else:                                   # codon G3: C34 and U34
            w[c] = ((1 - s[("C", "G")]) * tgcn.get(a_wc, 0)
                    + (1 - s[("U", "G")]) * tgcn.get(_revcomp(c[:2] + "A"), 0))
    if prokaryote:
        # AUA is decoded by tRNA-Ile2, whose C34 carries lysidine and reads A3.
        # Its anticodon CAU is the same string as tRNA-Met's, which is why the
        # parser keeps it under a separate key.
        w["ATA"] = w.get("ATA", 0.0) + (1 - s[("L", "A")]) * tgcn.get(ILE2_KEY, 0)
    return w


def tai(seq, trna_gene_copy_numbers, prokaryote=True, include_singletons=False,
        weights=None):
    """tRNA adaptation index (dos Reis, Savva & Wernisch 2004).

    Parameters
    ----------
    seq : str
        In-frame coding sequence.
    trna_gene_copy_numbers : dict
        Anticodon (uppercase DNA, for example ``"AGC"``) to gene copy number,
        as produced by :func:`trna_copy_numbers_from_tRNAscan`. The bacterial
        lysidine tRNA-Ile is expected under the key ``"CAT_Ile2"``.
    prokaryote : bool
        Apply the lysidine rule for AUA. Set ``False`` for eukaryotes.
    include_singletons : bool
        Include Met and Trp codons. dos Reis's implementation leaves them out,
        as CAI does, because they carry no synonymous choice; including them
        shifts the index by a constant factor that depends on the genome, so
        values computed with and without are not comparable.
    weights : dict, optional
        Precomputed relative adaptiveness, from a previous call. Supplying this
        avoids rebuilding the table for every gene in a genome.

    Returns
    -------
    float or None
        ``None`` when no tRNA data are available, when the table yields no
        non-zero weights, or when the sequence has no scoreable codons. A
        missing tAI is reported in the output table rather than imputed: the
        alternative is a number that looks like an answer and is not.

    Notes
    -----
    The absolute adaptiveness of codon i is
    ``W_i = sum_j (1 - s_ij) * tGCN_j`` over the tRNAs j that can decode it,
    with the selective constraints ``s`` of :data:`TAI_S_VALUES`. Weights are
    then normalised by the maximum, and codons with ``W_i = 0`` take the
    geometric mean of the non-zero weights, which is the correction in the
    original paper and keeps the geometric mean finite. tAI itself is the
    geometric mean of ``w`` over the codons of the gene.

    The ``s`` values were fitted to E. coli and S. cerevisiae expression data.
    They are used unchanged here, which is standard practice, but they are
    empirical constants from two organisms and not physical ones. Gene copy
    number is also only a proxy for tRNA abundance; it is a good proxy in fast
    growing bacteria and a poor one in slow growers and in eukaryotes.
    """
    if weights is None:
        if not trna_gene_copy_numbers:
            return None
        weights = tai_weights(trna_gene_copy_numbers, prokaryote=prokaryote)
    if not weights:
        return None
    total = n = 0
    for c in codon_list(seq):
        aa = CODON2AA.get(c)
        if aa is None or aa == "*":
            continue
        if not include_singletons and len(FAMILY[aa]) == 1:
            continue
        wi = weights.get(c)
        if not wi:
            continue
        total += math.log(wi)
        n += 1
    return math.exp(total / n) if n else None


def tai_weights(trna_gene_copy_numbers, prokaryote=True):
    """Normalised relative adaptiveness ``w`` for tAI, one value per sense codon."""
    if not trna_gene_copy_numbers:
        return {}
    tgcn = {str(k).upper().replace("U", "T"): float(v)
            for k, v in trna_gene_copy_numbers.items()}
    tgcn[ILE2_KEY] = float(trna_gene_copy_numbers.get(ILE2_KEY, 0))
    raw = _tai_weights(tgcn, prokaryote=prokaryote)
    mx = max(raw.values()) if raw else 0.0
    if mx <= 0:
        return {}
    w = {c: v / mx for c, v in raw.items()}
    nonzero = [v for v in w.values() if v > 0]
    if not nonzero:
        return {}
    geo = math.exp(sum(math.log(v) for v in nonzero) / len(nonzero))
    return {c: (v if v > 0 else geo) for c, v in w.items()}


_TRNA_AA = re.compile(r"^(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Ile2|Leu|Lys|"
                      r"Met|fMet|iMet|Phe|Pro|Ser|SeC|Sec|Thr|Trp|Tyr|Val|Pyl|"
                      r"Sup|Undet)$", re.I)
_ANTICODON = re.compile(r"^[ACGTUacgtu]{3}$")


def trna_copy_numbers_from_tRNAscan(path, include_pseudo=False,
                                    include_undetermined=False):
    """Count tRNA genes per anticodon from tRNAscan-SE output.

    Accepts the default tabular output of tRNAscan-SE 1.x and 2.x. The amino
    acid is taken from the ``Type`` column and the anticodon from the ``Anti
    Codon`` column, both located by shape rather than by fixed offset, because
    the column layout changed between versions.

    Anticodons are returned as uppercase DNA. Bacterial tRNA-Ile2, reported
    either as type ``Ile2`` or as ``Ile`` with anticodon ``CAT``, is stored
    under :data:`ILE2_KEY` rather than merged into the ``CAT`` count, because
    ``CAT`` has to mean tRNA-Met alone for the wobble rules to work.

    Suppressor and undetermined tRNAs are dropped by default, as are predicted
    pseudogenes, since neither contributes to decoding.
    """
    counts = defaultdict(float)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "-", "Sequence",
                                                    "Name", "--")):
                continue
            f = line.split()
            if len(f) < 6:
                continue
            aa = anti = None
            for i in range(len(f) - 1):
                if _TRNA_AA.match(f[i]) and _ANTICODON.match(f[i + 1]):
                    aa, anti = f[i], f[i + 1].upper().replace("U", "T")
                    break
            if aa is None:
                continue
            note = " ".join(f).lower()
            if not include_pseudo and "pseudo" in note:
                continue
            if not include_undetermined and aa.lower() in ("undet", "sup"):
                continue
            if aa.lower() == "ile2" or (aa.lower() == "ile" and anti == "CAT"):
                counts[ILE2_KEY] += 1
            else:
                counts[anti] += 1
    return dict(counts)


# ------------------------------------------------------- gene versus host
def _paired_frequencies(rscu_a, rscu_b):
    """Within-family frequencies for the families both tables can speak about.

    RSCU is a per-family quantity, so a table of RSCU values is not a
    probability distribution. Dividing by family size turns each family back
    into within-family frequencies summing to 1, and the families are then
    weighted equally. Equal weighting is a choice: weighting by amino-acid
    usage instead would make the distance partly a measure of protein
    composition, which is not what a codon-usage comparison is for.
    """
    pa, pb, n_families = [], [], 0
    for aa, fam in FAMILY.items():
        if len(fam) == 1:
            continue
        va = [rscu_a.get(c) for c in fam]
        vb = [rscu_b.get(c) for c in fam]
        if any(v is None for v in va) or any(v is None for v in vb):
            continue
        k = len(fam)
        pa.extend(v / k for v in va)
        pb.extend(v / k for v in vb)
        n_families += 1
    if not n_families:
        return [], [], 0
    pa = [v / n_families for v in pa]
    pb = [v / n_families for v in pb]
    return pa, pb, n_families


def cufs(rscu_a, rscu_b):
    """Codon usage frequency similarity, as a distance in [0, 1].

    Implements the Endres-Schindelin metric used by Diament, Pinter and Tuller
    (2014): the square root of the Jensen-Shannon divergence between two codon
    usage distributions. With base-2 logarithms the divergence is bounded by 1,
    so 0 means identical usage and 1 means the two never choose the same codon.

    Despite the name it is a distance, not a similarity. Use ``1 - cufs(...)``
    if a similarity is wanted, and say which was reported.

    Families empty in either table are skipped; the number of families actually
    compared is what limits how much a short gene can say. Returns ``None`` if
    no family is shared.

    .. warning::
       **This distance is dominated by gene length and must not be compared
       across genes of different length without a length control.** A short
       gene estimates each family's frequencies from few codons, so its RSCU
       table is noisy and its distance from any reference is inflated, with no
       difference in codon preference at all. Measured on the 4,016 ordinary
       chromosomal genes of *Escherichia coli* K-12 MG1655, where no mobility
       contrast exists by construction, Spearman's rho between CUFS to the
       genome centroid and coding length is -0.81 (n = 4,016). The same
       quantity for :func:`rscu_distance` is -0.77 for the euclidean metric and
       -0.69 for the correlation metric, against **+0.20** for CAI and -0.13
       for ENC.

       The sign of the CAI figure matters and an earlier version of this
       warning printed it as -0.20, which was wrong. CAI's weak length
       dependence runs the *other* way from the distances: a long gene tends
       to score slightly higher, not lower. So for a long gene the length
       artefact flatters its CAI and works against finding a CAI deficit,
       while it inflates its distance-based scores' apparent similarity to the
       host. Remeasured on the same 4,016 genes with independent code on
       2026-09-24, rho for CAI is +0.197 and for Fop +0.198 (see
       ``classes/01_colistin_mcr/scripts/210_rerun_housekeeping_contrast.py``).

       Part of the mechanism is countable: the number of amino-acid families a
       gene's RSCU table can be compared on at all is itself a length meter,
       rho +0.47 across the same 4,016 genes, running from 9 to 18 of the 18
       degenerate families. A short gene is scored on fewer families, and on
       each of them from fewer codons.

       The practical consequence is large: the 13 mcr colistin-resistance
       families, at 538 to 565 codons, score a *lower* raw CUFS to the E. coli
       centroid (0.202) than E. coli's own chromosomal resistance genes (0.249)
       purely because they are longer, which reverses the correct conclusion.
       Use a length-matched comparison, or an index with weak length dependence
       such as CAI or Fop, whenever lengths differ. Length-match without
       replacement: matching 13 mcr genes to a pool of 10 candidate partners
       with replacement leaves 7 distinct partners and turns a paired test's
       nominal 2**13 labellings into far fewer independent ones.
    """
    pa, pb, n = _paired_frequencies(rscu_a, rscu_b)
    if not n:
        return None
    jsd = 0.0
    for x, y in zip(pa, pb):
        m = 0.5 * (x + y)
        if m <= 0:
            continue
        if x > 0:
            jsd += 0.5 * x * math.log2(x / m)
        if y > 0:
            jsd += 0.5 * y * math.log2(y / m)
    return math.sqrt(max(jsd, 0.0))


def rscu_distance(rscu_a, rscu_b, metric="euclidean"):
    """Distance between two RSCU tables, over the families both define.

    Parameters
    ----------
    metric : {"euclidean", "manhattan", "cosine", "correlation"}
        ``"correlation"`` returns ``1 - r`` with Pearson ``r`` across codons,
        which ignores overall magnitude and is the least sensitive to a short
        gene having extreme RSCU values purely from small counts.

    Returns ``None`` if the tables share no usable family. A short gene leaves
    many families empty, so check how many codons went into the comparison
    before reading anything into the number.

    .. warning::
       Strongly length dependent, for the reason given under :func:`cufs`.
       Across the 4,016 ordinary chromosomal genes of *E. coli* K-12 MG1655,
       Spearman's rho against coding length is -0.77 for ``"euclidean"`` and
       -0.69 for ``"correlation"``. Never compare genes of different length on
       this scale without a length-matched control.
    """
    a, b = [], []
    for aa, fam in FAMILY.items():
        if len(fam) == 1:
            continue
        va = [rscu_a.get(c) for c in fam]
        vb = [rscu_b.get(c) for c in fam]
        if any(v is None for v in va) or any(v is None for v in vb):
            continue
        a.extend(va)
        b.extend(vb)
    if not a:
        return None
    if metric == "euclidean":
        return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))
    if metric == "manhattan":
        return sum(abs(x - y) for x, y in zip(a, b))
    if metric == "cosine":
        na = math.sqrt(sum(x * x for x in a))
        nb = math.sqrt(sum(y * y for y in b))
        if na == 0 or nb == 0:
            return None
        return 1.0 - sum(x * y for x, y in zip(a, b)) / (na * nb)
    if metric == "correlation":
        n = len(a)
        ma, mb = sum(a) / n, sum(b) / n
        sa = math.sqrt(sum((x - ma) ** 2 for x in a))
        sb = math.sqrt(sum((y - mb) ** 2 for y in b))
        if sa == 0 or sb == 0:
            return None
        r = sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb)
        return 1.0 - r
    raise ValueError("unknown metric %r" % (metric,))
