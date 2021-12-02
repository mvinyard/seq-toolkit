
# _SequenceGenerator.py

__module_name__ = "_SequenceGenerator.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import numpy as np

def _set_weight_simplex(A=1, C=1, G=1, T=1):
    
    """
    Change the composition of weights for sampling bases at random.

    Parameters:
    ----------
    N {A, C, G, T}
        proportions of bases to be sampled. simplex. 

    Returns:
    -------
    Initializes self.Gene and creates/modifies self.SeqGen

    Notes:
    ------
    (1) As an example, if A=2, the simplex vector would appear as: [0.4, 0.2, 0.2, 0.2]
    """
    
    bases = np.array([A, C, G, T])
    return bases / bases.sum()

class _SequenceGenerator:
    
    def __init__(self, A=1, C=1, G=1, T=1,):

        """
        Initialize random sequence generator.

        Parameters:
        ----------
        N {A, C, G, T}
            proportions of bases to be sampled. simplex. 

        Returns:
        -------
        self.bases
            Array of available bases.
            type: numpy.ndarray

        self.weights
            Weights to be passed for base selection.
            type: numpy.ndarray

        Notes:
        ------
        (1) As an example, if A=2, the simplex vector would appear as: [0.4, 0.2, 0.2, 0.2]
        """
        
        self.bases = np.array(["A", "C", "G", "T"])
        self.weights = _set_weight_simplex(A, C, G, T)      
        
    def simulate(self, n_bases, return_seq=True):
        
        
        """"
        Simulate a DNA sequence of arbitrary length.
        
        Parameters:
        -----------
        n_bases
            Number of bases to be selected.
            
        return_seq
            Indicates if the generated sequence should be returned directly. Saved as a subclass object by default.
            type: bool
            default: True
        
        Returns:
        --------
        [ optional ] self.seq
            Generated sequence of length(n_bases)
            type: str
        
        Notes:
        ------
        """
          
        self.seq = "".join(np.random.choice(self.bases, n_bases, p=self.weights))
        
        if return_seq:
            return self.seq
