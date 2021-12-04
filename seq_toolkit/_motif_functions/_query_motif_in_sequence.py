
# _query_motif_in_sequence.py

__module_name__ = "_query_motif_in_sequence.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import licorice
import pandas as pd
import regex


# local imports #
# ------------- #
from ._isolate_constraining_sequence_motif import _isolate_constraining_sequence_motif
from .._sequence_functions._SequenceManipulation import _SequenceManipulation


def _query_motif_in_sequence(sequence, motif, strand, start_key, end_key, verbose):
    
    """"""
    
    if verbose:
        strand_str = licorice.font_format(strand, ["BOLD"])
        motif_str = licorice.font_format(motif, ["BOLD", "GREEN"])
        print("Searching the {} strand of the provided sequence for: {} ...".format(strand_str, motif_str))
    
    Motifs = {}
    motif_matches = regex.finditer(motif, sequence, overlapped=True)
    
    for n, match in enumerate(motif_matches):
        motif_position = match.span()
        Motifs[n] = {}
        Motifs[n][start_key] = motif_position[0] - 1
        Motifs[n][end_key]   = motif_position[1]
        
    motif_df = pd.DataFrame(Motifs).T
    
    motif_key = start_key.split(".")[0]
    motif_df["{}.strand".format(motif_key)] = strand
    
    return motif_df


def _query_motif_bistrand(sequence, motif, motif_key=False, verbose=True):
    
    """
    Look for a motif in both strands of a given DNA sequence. 
    
    Parameters:
    -----------
    sequence
        String to be searched for a motif.
        type: str
    
    motif
        String to be searched as a sub-string of the provided `sequence`
        type: str
        
    motif_key
        String to indicate the colunmn title for the motif start and end.
        type: str
    
    verbose
        Indicates if messages to the user should be printed or silenced. 
        type: bool
        default: True
    
    Returns:
    --------
    motif_df
        Pandas DataFrame containing all occurances of the searched motif.
        type: pandas.DataFrame
        
    
    Notes:
    ------
    """
    
    if not motif_key:
        motif_key="motif"
    
    
    start_key="{}.start".format(motif_key)
    end_key="{}.end".format(motif_key)
    
    searchable_motif = _isolate_constraining_sequence_motif(motif, verbose=verbose)
    
    Sequence = _SequenceManipulation(sequence)
    Sequence.reverse_complement()        

    pos_df = _query_motif_in_sequence(sequence, 
                                      searchable_motif, 
                                      "+", 
                                      start_key, 
                                      end_key,
                                      verbose,)
    
    neg_df = _query_motif_in_sequence(Sequence.reverse_complement_sequence, 
                                      searchable_motif, 
                                      "-", 
                                      start_key, 
                                      end_key,
                                      verbose,)        

    neg_df[start_key] = len(sequence) - neg_df[start_key]
    neg_df[end_key] = len(sequence) - neg_df[end_key]
    
    motif_df = pd.concat([pos_df, neg_df]).sort_values([start_key]).reset_index(drop=True)
    
    
    
    if verbose:
        print("\nIdentified {} instances of the {} motif in the provided sequence.".format(len(motif_df), searchable_motif))
    
    return motif_df