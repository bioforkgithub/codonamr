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
#: three silently rejects real genes: measured on 2026-09-24, the 3-codon set
#: rejects 18 of the 9,786 CDS in shared/data/AMR_CDS.fa and 0 to 47 CDS per
#: genome across the five host genomes in classes/01_colistin_mcr/data/
#: host_cache, all of them with a legal table-11 initiator.
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
    """Split an in-frame DNA string into codons."""
    return [seq[i:i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
