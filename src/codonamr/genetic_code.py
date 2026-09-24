"""Standard genetic code tables and synonymous families.

Codon order follows the TCAG convention used by CodonW and by most codon-usage
literature, so tables printed by this package line up with published ones.
"""
from collections import defaultdict

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON2AA = dict(zip(CODONS, AAS))
STOPS = frozenset(c for c, a in CODON2AA.items() if a == "*")

#: Initiation codons of NCBI translation table 11 (Bacterial, Archaeal and
#: Plant Plastid), which is the code this package is for. Table 11 marks seven
#: codons as initiators: ATG, GTG, TTG, ATT, ATC, ATA and CTG. All seven are
#: translated as Met when they occur in the initiator position, which is why
#: AMRFinderPlus and Prokka call a protein-level 100 per cent allele match on a
#: CDS whose first codon is, for example, ATA. Restricting QC to the first
#: three rejects real genes. Measured on 2026-09-24 (reproduce with
#: classes/01_colistin_mcr/scripts/141_referee_check_cds_start_set.py), the
#: 3-codon set wrongly rejects 18 of the 9,786 CDS in shared/data/AMR_CDS.fa,
#: every one of which carries a legal table-11 initiator, and 0 to 32 CDS per
#: genome across the five host genomes in classes/01_colistin_mcr/data/
#: host_cache (Enterobacter hormaechei 32, Acinetobacter baumannii 13,
#: Salmonella Typhimurium LT2 12, E. coli O157:H7 Sakai 2, Klebsiella
#: pneumoniae HS11286 0). Those per-genome figures are the CDS recovered by
#: the wider set, not the bad_start totals: the two draft assemblies also
#: carry CDS that begin at no initiator at all (15 in E. hormaechei, 9 in
#: A. baumannii) and those stay rejected, correctly.
STARTS_TABLE11 = ("ATG", "GTG", "TTG", "ATT", "ATC", "ATA", "CTG")

#: The conservative 3-codon subset used before 2026-09-24. Kept so that a
#: caller can reproduce the older QC counts exactly.
STARTS_CONSERVATIVE = ("ATG", "GTG", "TTG")

_fam = defaultdict(list)
for _c, _a in CODON2AA.items():
    if _a != "*":
        _fam[_a].append(_c)
FAMILY = dict(_fam)

#: the 18 amino acids with more than one codon (excludes Met, Trp and stops)
DEGENERATE = {a: cs for a, cs in FAMILY.items() if len(cs) > 1}

#: amino acids encoded by a single codon, which carry no synonymous information
SINGLETON = {a for a, cs in FAMILY.items() if len(cs) == 1}


def codon_list(seq):
    """Split an in-frame DNA string into codons.

    The sequence is upper-cased first. Without that, a lower-case FASTA (which
    ``shared/data/AMR_CDS.fa`` is) silently produced ``None`` from ``gc3s``,
    ``enc``, ``scuo``, ``fop`` and ``cai`` and all-zero counts from
    ``codon_counts``, while ``gc_content``, which upper-cases on its own, kept
    returning a correct value. That mixture is worse than an error. Callers
    that go through ``codonamr.qc.check_cds`` were never affected, since it
    upper-cases and returns the cleaned sequence.
    """
    seq = seq.upper()
    return [seq[i:i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
