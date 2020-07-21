#!/bin/bash

import os
import pandas as pd
import numpy as np


homepath = os.path.expanduser("~")


def concatenate_files(*file_paths):
    """
    :param file_paths: Paths to files to concatenate.
    :return: pandas.DataFrame of estimated counts per sample.
    """
    for file_number, file_path in enumerate(file_paths, 1):
        filename = os.path.basename(os.path.dirname(file_path))
        new_or_old = os.path.basename(os.path.dirname(os.path.dirname(file_path))).lower()

        if file_number == 1:
            abundance_file = pd.read_csv(file_path, sep="\t")[["target_id", "est_counts"]]
            abundance_file.columns = ["target_id", f"est_counts_{filename}_{new_or_old}"]
        else:
            abundance_file = pd.read_csv(file_path, sep="\t")["est_counts"]
            abundance_file = abundance_file.rename(f"est_counts_{filename}_{new_or_old}")
        
        try:
            counts = pd.concat([counts, abundance_file], axis=1)
        except UnboundLocalError:
            counts = abundance_file

    # Replacing NaN with zeros and converting floats to integers for following DESeq2 analysis
    counts.iloc[:, list(range(1, len(file_paths) + 1))] = counts.iloc[:, list(range(1, len(file_paths) + 1))].fillna(0).astype(int)
    return counts


if __name__ == "__main__":
    cancer_samples = concatenate_files(
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/New/ERR358487/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/New/ERR358488/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/Old/ERR358487/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/Old/ERR358488/abundance.tsv",
    )
    cancer_samples.to_csv(f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Files/cancer_kallisto_counts.csv", index=False)

    normal_samples = concatenate_files(
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/New/ERR358485/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/New/ERR358486/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/Old/ERR358485/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/Old/ERR358486/abundance.tsv",
    )
    normal_samples.to_csv(f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Files/normal_kallisto_counts.csv", index=False)

