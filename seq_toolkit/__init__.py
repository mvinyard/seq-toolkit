# main __init__.py

__module_name__ = "__init__.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


from ._gene_functions._GeneGenerator import _GeneGenerator as Gene

from ._motif_functions._query_motif_in_sequence import _query_motif_bistrand as query_motif
from ._motif_functions._isolate_constraining_sequence_motif import _isolate_constraining_sequence_motif as isolate_searchable_motif

from ._sequence_functions._SequenceManipulation import _SequenceManipulation as SequenceManipulator
from ._sequence_functions._SequenceGenerator import _SequenceGenerator as Seq

from ._fetch_chromosome_sequence import _fetch_chromosome_sequence as fetch_chromosome
