
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import vinplots


from ._SequenceGenerator import _SequenceGenerator

def _define_gene_exons(n_bases, n_boundaries, boundary_spacing):
    
    """
    Choose positions for gene exons. 
    
    Parameters:
    -----------
    n_bases
        Length of matching inputs sequence
        type: int
        
    n_boundaries
        Number of samples drawn from the range of len(sequence)
        type: int
    
    boundary_spacing
        choice of every "n" boundary
        type: int
    
    Returns:
    --------
    exon_df
    
    Notes:
    ------
    """
    
    ExonDict = {}
    ExonDict['Start'] = []
    ExonDict['End'] = []
    
    exon_bounds = np.sort(np.random.choice(n_bases, n_boundaries))
    
    for i in range(len(exon_bounds)-1):
        if i % boundary_spacing == 0:
            ExonDict['Start'].append(exon_bounds[i])
            ExonDict['End'].append(exon_bounds[i + 1])
           
    exon_df = pd.DataFrame.from_dict(ExonDict)
    
    return exon_df

def _construct_gene_plot(gene_width, n_ticks):
    
    """
    Create the framework for the gene figure. 
    
    Parameters:
    -----------
    gene_width
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
    xt = ax.set_xticks(np.linspace(0, gene_width, n_ticks))
    yt = ax.set_yticks([])
    
    return fig, ax

def _plot_gene(exon_df, color="navy", save=False):
    
    """
    Plot a simulated gene. 

    Parameters:
    -----------
    color
        type: str
        default: "navy"
        
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
    """
    
    gene_width = exon_df.End.max() - exon_df.Start.min()
    
    fig, ax = _construct_gene_plot(gene_width)
    
    plt.hlines(1, exon_df.Start.min(), exon_df.End.max(), color=color, zorder=2)
    plt.ylim(.95, 1.1)
    
    for i, exon in exon_df.iterrows():
        plt.vlines(exon.Start, 0.95, 1.025, color="lightgrey", linestyle="--", lw=1)
        plt.vlines(exon.End, 0.95, 1.025, color="lightgrey", linestyle="--", lw=1)
        plt.text(x=exon.Start + gene_width/500, y=1.03, s="{}-{}".format(exon.Start, exon.End), rotation=25)
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
        
        
    def create(self, n_bases, n_boundaries=30, boundary_spacing=5, return_gene=True):
        
        self.n_bases = n_bases
        self.Gene["seq"] = self.seq = self.SeqGen.simulate(n_bases, return_seq=True)
        self.Gene["exons"] = self.exon_df = _define_gene_exons(n_bases, n_boundaries, boundary_spacing)
        
        
        if return_gene:
            return self.seq
        
    def plot(self, color="navy", save=False):
        
        """
        Plot a simulated gene. 
        
        Parameters:
        -----------
        color
            type: str
            default: "navy"
            
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
        
        _plot_gene(self.exon_df, color, save)
