# _merge_reduce_features.py

__module_name__ = "_merge_reduce_features.py"
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

        chrom_vals = self._grouped_df["Chromosome"]
        start_vals = self._grouped_df["Start"].aggregate(np.min)
        end_vals = self._grouped_df["End"].aggregate(np.max)

        self.merged_df = (
            pd.DataFrame([start_vals, end_vals])
            .T.reset_index()
            .merge(self._clustered_df, on=["Cluster", "Start", "End"])
        ).drop("Cluster", axis=1)[["Chromosome", "Start", "End"]]

    def write_bed(self, out_path="merged_features.bed"):

        """

        Can be downloaded and visualized directly in IGV.

        """

        self.merged_df.to_csv(out_path, sep="\t", header=False, index=False)
