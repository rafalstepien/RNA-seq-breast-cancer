#!/bin/bash
import os
import pandas as pd
import numpy as np

def get_difference(dataframe):
    dataframe['diff_new'] = abs(dataframe[:, 1] - dataframe[:, 2])
    dataframe['diff_old'] = abs(dataframe[:, 3] - dataframe[:, 4])
    return dataframe

def concatenate_files(*file_paths):
    """
    :param file_paths: Paths to files to concatenate.
    :return: pandas.DataFrame of estimated counts per sample.
    """
    for _, file_path in enumerate(file_paths, 1):

        # Get name of a file from path
        filename = os.path.basename(os.path.dirname(file_path))

        # Get name of a directory above (Old/New)
        new_or_old = os.path.basename(os.path.dirname(os.path.dirname(file_path))).lower()

        # Read file, rename columns and drop transcript version
        abundance_file = pd.read_csv(file_path, sep="\t")[["target_id", "est_counts"]]
        abundance_file.columns = ["target_id", f"est_counts_{filename}_{new_or_old}"]
        abundance_file.target_id = pd.Series([row.split('.')[0] for _, row in abundance_file.target_id.items()])

        try:
            counts = counts.merge(abundance_file, on='target_id')
        except UnboundLocalError:
            counts = abundance_file

    # Replacing NaN with zeros and converting floats to integers for following DESeq2 analysis
    counts.iloc[:, list(range(1, len(file_paths) + 1))] = counts.iloc[:, list(range(1, len(file_paths) + 1))].fillna(0).astype(int)
    # counts = get_difference(counts)
    return counts


if __name__ == "__main__":
    homepath = os.path.expanduser("~")
    cancer_samples = concatenate_files(
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/New/ERR358487/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/New/ERR358488/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/Old/ERR358487/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/Old/ERR358488/abundance.tsv",
    )
    cancer_samples.to_csv(f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Outputs/concatenate_files_output/cancer_kallisto_counts.csv", index=False)

    normal_samples = concatenate_files(
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/New/ERR358485/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/New/ERR358486/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/Old/ERR358485/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/Old/ERR358486/abundance.tsv",
    )
    normal_samples.to_csv(f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Outputs/concatenate_files_output/normal_kallisto_counts.csv", index=False)
    