
__module_name__ = "_GenomicFeatures.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(
    [
        "vinyard@g.harvard.edu",
    ]
)


# package imports #
# --------------- #
import numpy as np
import pandas as pd
import os

try:
    import pyranges
except:
    os.system("sudo apt-get install gcc")
    os.system("pip install pyranges")
    import pyranges


def _cluster_df_features(df):

    """
    Requires the standard notation: df[['Chromosome', 'Start', 'End']]
    """

    ClusteredRanges_df = (
        pyranges.PyRanges(chromosomes=df.Chromosome, starts=df.Start, ends=df.End)
        .cluster()
        .as_df()
        .drop_duplicates()
    )

    return ClusteredRanges_df


def _cluster_overlapping_features(df):

    """"""

    clustered_df = _cluster_df_features(df)
    grouped_df = clustered_df.groupby("Cluster")

    return clustered_df, grouped_df


class _GenomicFeatures:

    """general module for merge-reducing a pandas DataFrame with start and stop feature designations."""

    def __init__(self, df):

        """"""

        self.df = df
        self._clustered_df, self._grouped_df = _cluster_overlapping_features(df)

    def merge(self):

        self._chrom_vals = self._grouped_df["Chromosome"].aggregate(np.unique).values
        self._start_vals = self._grouped_df["Start"].min()
        self._end_vals = self._grouped_df["End"].max()

        self.merged_df = pd.DataFrame(
            [self._chrom_vals, self._start_vals, self._end_vals],
            index=["Chromosome", "Start", "End"],
        ).T

    def write_bed(self, out_path="merged_features.bed"):

        """
        Can be downloaded and visualized directly in IGV.
        """

        self.merged_df.to_csv(out_path, sep="\t", header=False, index=False)