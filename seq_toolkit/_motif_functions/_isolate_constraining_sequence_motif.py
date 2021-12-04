
# _isolate_constraining_sequence_motif.py

__module_name__ = "_isolate_constraining_sequence_motif.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import licorice


def _isolate_constraining_sequence_motif(motif, verbose=False):
    
    """
    Get rid of the Ns in the motif to be queried.
    
    Parameters:
    -----------
    motif
        Motif to be stripped of any "N" values.
        type: str
    
    verbose
        type: bool
        default: False
    
    Notes:
    ------
    (1) Keeps N values that are in-between other bases. In other
        words, only N values that are on the end of sequences are
        trimmed.
    
    (2) To be implemented: Searchable motif generatator that subs 
        {A, G, T, C} for N wherein N is searchable / in the middle 
        of the motif. 
    """
    
    motif_ = motif.strip('N')
    
    if verbose:
        m = licorice.font_format(motif, ['BOLD'])
        m_=licorice.font_format(motif_, ['BOLD', 'GREEN'])
        print("\nSearching for motif: {} (from {})\n".format(m_, m,))
    
    return motif_