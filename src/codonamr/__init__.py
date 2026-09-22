"""codonamr: codon usage analysis of antimicrobial resistance genes.

Supply sequence identifiers; the package fetches the coding sequences, works
out which genome they came from, builds that host's codon table from its RefSeq
reference genome, computes the standard codon-usage indices, and then applies
the controls that decide whether any difference it finds is real: a
length-matched null, because the indices are biased on short genes, and a
partial correlation on per-nucleotide composition, because base composition
drives codon choice on its own. The controls are themselves validated by
injecting an association of known size and by injecting none.

Quickstart
----------
>>> from codonamr import analyse
>>> result = analyse(ids=["NG_048025.1"], out_dir="results",
...                  email="you@example.org")        # doctest: +SKIP

or from a shell, with CODONAMR_EMAIL set::

    codonamr run --ids ids.txt --out results/

Primary references for the implemented methods are named in each module and
collected in docs/metrics.md.
"""
__version__ = "0.1.0"

from .controls import (base_fractions, composition_control,
                       false_positive_rate, injection_recovery,
                       length_matched_null, spearman_partial)
from .fetch import CdsRecord, fetch_cds, read_fasta, resolve_host_genome
from .folding import available as folding_available
from .folding import dg_windowed
from .genetic_code import CODON2AA, CODONS, DEGENERATE, FAMILY, codon_list
from .host import (cufs, host_codon_table, optimal_codons, reference_set,
                   rscu_distance, tai, trna_copy_numbers_from_tRNAscan)
from .metrics import (aa_composition, build_w, cai, codon_counts, enc, fop,
                      gc3s, gc_content, rscu, scuo)
from .pipeline import analyse
from .qc import REASONS, check_cds, clean_cds
from .randomise import (junction_counts, randomise_replace, randomise_shuffle,
                        randomise_shuffle_dinuc)

__all__ = [
    "__version__", "analyse",
    # sequence input
    "CdsRecord", "fetch_cds", "read_fasta", "resolve_host_genome",
    # genetic code
    "CODONS", "CODON2AA", "FAMILY", "DEGENERATE", "codon_list",
    # quality control
    "check_cds", "clean_cds", "REASONS",
    # indices
    "codon_counts", "aa_composition", "rscu", "enc", "build_w", "cai",
    "gc_content", "gc3s", "fop", "scuo",
    # host comparison
    "host_codon_table", "reference_set", "optimal_codons", "tai",
    "trna_copy_numbers_from_tRNAscan", "cufs", "rscu_distance",
    # nulls
    "randomise_replace", "randomise_shuffle", "randomise_shuffle_dinuc",
    "junction_counts",
    # controls and their validation
    "length_matched_null", "composition_control", "spearman_partial",
    "base_fractions", "injection_recovery", "false_positive_rate",
    # folding
    "dg_windowed", "folding_available",
]
