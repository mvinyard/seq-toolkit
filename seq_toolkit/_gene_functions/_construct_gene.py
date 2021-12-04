
# _construct_gene.py

__module_name__ = "_construct_gene.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import pandas as pd
import numpy as np

def _construct_repetitive_feature(feature_space_sum, 
                                  n_features, 
                                  feature, 
                                  start_key, 
                                  end_key,
                                  feature_key,
                                  zero_start=False,):
    
    """
    Given a length of space to occupy, a set number of features is generated. 
    
    Parameters:
    -----------
    feature_space_sum
        Total space from which feature boundaries are selected.
        type: int
    
    n_features
        Number of feature boundaries constructed within the feature_space_sum
        type: int
    
    feature
        Name of the feature to be denoted.
        type: str
    
    zero_start
        Indicates whether the first feature should start at 0 or not.
        type: bool
        default: False
    
    Returns:
    --------
    feature_df
    
    Notes:
    ------
    (1) Example: the following parameters:
            feature_space_sum = 500
            n_features = 20
            zero_start = False
        
        Will return a DataFrame wherein each row is a feature (there will be 20 features/rows)
        and the boundaries are denoted by randomly choosing integers betweeen 0 and 500. Since
        `zero_start` is set to `False`, the first feature boundary may not necessarily = zero.
    """
    
    FeatureDict = {}
    FeatureDict[start_key] = []
    FeatureDict[end_key]   = []
    FeatureDict[feature_key] = feature    
    feature_dividers = np.sort(np.random.randint(0, feature_space_sum, n_features))
    
    if zero_start:
        feature_dividers[0] = 0

    for i in range(len(feature_dividers)):
        FeatureDict[start_key].append(feature_dividers[i])
        if i == len(feature_dividers)-1:
            FeatureDict[end_key].append(feature_space_sum)
        else:
            FeatureDict[end_key].append(feature_dividers[i+1])
    
    feature_df = pd.DataFrame(data=FeatureDict)
    
    return feature_df

def _preassemble_genebody_df(intron_df, exon_df):

    rows = []
    rows.append(exon_df.iloc[0])
    for i in range(1, len(intron_df) + 1):
        try:
            rows.append(intron_df.iloc[i-1])
            rows.append(exon_df.iloc[i])
        except:
            rows.append(exon_df.iloc[i])
    
    pregene_df = pd.DataFrame(rows)
    
    return pregene_df

def _subtract_initial_UTR_space(gene_df, start_key, end_key):
    
    """"""
    
    inital_space = gene_df.iloc[0][start_key]
    
    gene_df[start_key] = gene_df[start_key] - inital_space
    gene_df[end_key] = gene_df[end_key] - inital_space
    
    return gene_df, inital_space

def _annotate_UTRs(gene_df, gene_length, feature_key, start_key, end_key):
    
    """
    
    Parameters:
    -----------
    
    Returns:
    --------
    
    """
    
    if gene_df.iloc[0][start_key] != 0:
        gene_df, total_UTR_space = _subtract_initial_UTR_space(gene_df, start_key, end_key)
        
        UTR = {}
        UTR["5prime"] = {}
        UTR["3prime"] = {}
        
        UTR_5p_len = int(total_UTR_space/np.random.randint(2,6))
        UTR_3p_len = total_UTR_space - UTR_5p_len
        
        # now add back the adjusted 5' UTR 
        gene_df[start_key] = gene_df[start_key] + UTR_5p_len
        gene_df[end_key] = gene_df[end_key] + UTR_5p_len
        
        UTR["5prime"][start_key] = 0
        UTR["5prime"][end_key] = UTR_5p_len
        UTR["3prime"][start_key] = gene_df.iloc[-1][end_key]
        UTR["3prime"][end_key] = gene_length
        
        UTR_df = pd.DataFrame(UTR).T
        UTR_df[feature_key] = "UTR"
        
        gene_df = pd.concat([gene_df, UTR_df]).sort_values([start_key]).reset_index(drop=True)
        return gene_df
        
    else:
        return gene_df
       

def _assemble_gene_df(intron_df, exon_df, start_key, end_key, feature_key):
    
    """"""
    
    pregene_df = _preassemble_genebody_df(intron_df, exon_df)
    
    contextualized_features = {}
    contextualized_features[start_key] = []
    contextualized_features[end_key]= []
    contextualized_features[feature_key] = []
    for i in range(len(pregene_df)):
        if i != 0:
            contextualized_features[start_key].append(pregene_df.iloc[i][start_key] + pregene_df.iloc[i-1][end_key],)
            contextualized_features[end_key].append(pregene_df.iloc[i][end_key] + pregene_df.iloc[i-1][end_key])
            contextualized_features[feature_key].append(pregene_df.iloc[i][feature_key])
        else:
            contextualized_features[start_key].append(pregene_df.iloc[i][start_key])
            contextualized_features[end_key].append(pregene_df.iloc[i][end_key])
            contextualized_features[feature_key].append(pregene_df.iloc[i][feature_key])
    
    gene_df = pd.DataFrame(contextualized_features).reset_index(drop=True)
    
    return gene_df

def _construct_gene(gene_length=50000, 
                    n_exons=15, 
                    min_exon_length=50, 
                    max_exon_length=2500, 
                    start_key="gene_feature.start",
                    end_key="gene_feature.end",
                    feature_key="gene_feature",
                    zero_start=False, 
                    verbose=False):
    
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
    gene_df
    
    Notes:
    ------
    """
    
    n_introns = n_exons-1
    exon_lengths = np.random.randint(min_exon_length, max_exon_length, n_exons)
    exon_sum = exon_lengths.sum()
    intron_sum = gene_length - exon_sum
    
    if verbose:
        print("Total exon length:\t{}\nTotal intron length:\t{}".format(exon_sum, intron_sum))
    
    exon_df = _construct_repetitive_feature(exon_sum, 
                                            n_exons, 
                                            "exon", 
                                            start_key, 
                                            end_key,
                                            feature_key,
                                            zero_start)
    
    intron_df =_construct_repetitive_feature(intron_sum, 
                                             n_introns, 
                                             "intron", 
                                             start_key, 
                                             end_key,
                                             feature_key,
                                             zero_start,)

    gene_df = _assemble_gene_df(intron_df, exon_df, start_key, end_key, feature_key)
    gene_df = _annotate_UTRs(gene_df, gene_length, feature_key, start_key, end_key)
    gene_df['length'] = gene_df[end_key] - gene_df[start_key]

    return exon_df, intron_df, gene_df