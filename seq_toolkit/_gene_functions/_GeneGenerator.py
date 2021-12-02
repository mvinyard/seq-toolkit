
# _GeneGenerator.py

__module_name__ = "_GeneGenerator.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import vinplots


# local imports #
# ------------- #
from .._sequence_functions._SequenceGenerator import _SequenceGenerator
from ._construct_gene import _construct_gene


def _define_gene_exons(gene_length, n_exons=15, min_exon_length=50, max_exon_length=2500):
    
    """
    Choose positions for gene exons. 
    
    Parameters:
    -----------
    gene_length
        Length of input gene body.
        type: int
        
    n_exons
        Number of desired exons.
        type: int
        default: 15
        
    min_exon_length
        Minimum exon length
        type: int
        default: 50
        
    max_exon_length
        Maximum exon length
        type: int
        default: 2500
    
    Returns:
    --------
    exon_df
    
    Notes:
    ------
    """

    exon_lengths = np.random.randint(min_exon_length, max_exon_length, n_exons)
    exon_start_positions = np.sort(np.random.randint(0, gene_length-exon_lengths[-1], n_exons))
    exon_end_positions = exon_start_positions + exon_lengths

    exon_df = pd.DataFrame(data = {'Start': exon_start_positions, 'End': exon_end_positions})

    return exon_df

def _construct_gene_plot(n_bases, n_ticks):
    
    """
    Create the framework for the gene figure. 
    
    Parameters:
    -----------
    n_bases
        Length of the gene being plotted. 
    
    n_ticks
        Number of ticks along the x-axis. 
        type: int
        
    Returns:
    --------
    fig, ax
    """

    fig = vinplots.Plot()
    fig.construct(nplots=1, ncols=1, figsize_width=2.5)
    fig.modify_spines(ax="all", spines_to_delete=['top', 'right', 'left'])
    ax = fig.AxesDict[0][0]
    xt = ax.set_xticks(np.linspace(0, n_bases, n_ticks))
    yt = ax.set_yticks([])
    
    return fig, ax

def _plot_gene(n_bases, gene_df, color="navy", n_ticks=11, save=False):
    
    """
    Plot a simulated gene. 

    Parameters:
    -----------
    n_bases
        Length of plotted gene.
        
    gene_df
        pandas DataFrame with columns: ['Start', 'End'] denoting the boundaries of the gene UTRs, introns, and exons.
        type: pandas.DataFrame
        
    color
        type: str
        default: "navy"
        
    n_ticks
        Number of ticks along the gene plot's x-axis.
        type: int
        default: 11
        
    save
        If not False, pass a string, which will trigger the object to save with figname=`save`.
        type: bool or str
        default: False
        
    Returns:
    --------
    None, prints (and/or saves) plot.

    Notes:
    ------
    (1) requires prior running of `Gene.create()`
    (2) 500 is conveniently/arbitrarily chosen as a scalar for offsetting the text-labels of exon spans. 
        May not suit all use-cases and could be updated in the future if needed. 
    """    
    
    fig, ax = _construct_gene_plot(n_bases, n_ticks)
    plt.ylim(.95, 1.1)
    
    plt.hlines(1, gene_df.Start.min(), gene_df.End.max(), color=color, zorder=2)
    text_offset = (gene_df.End.max() - gene_df.Start.min()) / 500
    
    exon_df = gene_df.loc[gene_df.Feature == "exon"]
    
    for i, exon in exon_df.iterrows():
        plt.vlines(exon.Start, 0.95, 1.025, color="lightgrey", linestyle="--", lw=1)
        plt.vlines(exon.End, 0.95, 1.025, color="lightgrey", linestyle="--", lw=1)
        plt.text(x=exon.Start + text_offset, y=1.03, s="{}-{}".format(exon.Start, exon.End), rotation=25)
        plt.hlines(1, exon.Start, exon.End, color=color, lw=15, zorder=2)
    
    if save:
        fig.savefig(save)

class _GeneGenerator:
    
    def __init__(self, A=1, C=1, G=1, T=1):
        
        """
        Initializes the `Gene` class.
        
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
        
        self.Gene = {}
        self.SeqGen = _SequenceGenerator(A, C, G, T)
        
        
    def create(self, 
               gene_length=50000,
               n_exons=15, 
               min_exon_length=50, 
               max_exon_length=2500, 
               return_gene=True, 
               zero_start=False,
               verbose=False,):
        
        """
        Executes creation of the gene within the set parameters.
        
        Parameters:
        -----------
        gene_length
            Length of input gene body.
            type: int

        n_exons
            Number of desired exons.
            type: int
            default: 15

        min_exon_length
            Minimum exon length
            type: int
            default: 50

        max_exon_length
            Maximum exon length
            type: int
            default: 2500

        return_gene
            Indicates if the gene sequence (self.seq) should be returned directly.
            By default, the gene sequence is also stored in memory as the subclass, self.seq.
            type: bool
            default: False
        
        
        Returns:
        --------
        [ optional ] self.seq
            Gene sequence
            type: str
        
        self.n_bases
        
        self.Gene["seq"]
            Gene sequence of len(n_bases)
            type: str
            
        self.Gene["exons"], self.exon_df
            pandas DataFrame with columns: ['Start', 'End'] denoting the boundaries of the gene exons.
            type: pandas.DataFrame
        
        Notes:
        ------
        """
        
        self.gene_length = gene_length
        self.Gene["seq"] = self.seq = self.SeqGen.simulate(gene_length, return_seq=True)
        self.gene_df = _construct_gene(gene_length, 
                                          n_exons, 
                                          min_exon_length, 
                                          max_exon_length, 
                                          zero_start, 
                                          verbose)
        
        
        if return_gene:
            return self.seq
        
    def plot(self, color="navy", n_ticks=11, save=False):
        
        """
        Plot a simulated gene. 
        
        Parameters:
        -----------
        color
            type: str
            default: "navy"
            
        n_ticks
            Number of ticks along the gene plot's x-axis.
            type: int
            default: 11
        
        save
            If not False, pass a string, which will trigger the object to save with figname=`save`.
            type: bool or str
            default: False
        
        Returns:
        --------
        None, prints plot.
        
        Notes:
        ------
        (1) `save` is not yet implemented
        (2) requires prior running of `Gene.create()`
        """
        
        _plot_gene(self.gene_length, self.gene_df, color, n_ticks, save)
