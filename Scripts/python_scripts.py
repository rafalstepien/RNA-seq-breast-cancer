#!/bin/bash

import pandas as pd
import os


def concatenate_files(*file_paths):
    """
    :param file_paths: Paths to files to concatenate.
    :return: pandas.DataFrame of estimated counts per sample.
    """
    for file_number, file_path in enumerate(file_paths, 1):
        filename = os.path.basename(os.path.dirname(file_path))

        if file_number == 1:
            abundance_file = pd.read_csv(file_path, sep="\t")[["target_id", "est_counts"]]
            abundance_file.columns = ["target_id", f"est_counts_{filename}"]
        else:
            abundance_file = pd.read_csv(file_path, sep="\t")["est_counts"]
            abundance_file = abundance_file.rename(f"est_counts_{filename}")
        
        try:
            counts = pd.concat([counts, abundance_file], axis=1)
        except UnboundLocalError:
            counts = abundance_file

    return counts


if __name__ == "__main__":
    homepath = os.path.expanduser("~")

    cancer_samples = concatenate_files(
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/New/ERR358487/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/New/ERR358488/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/Old/ERR358487/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/Old/ERR358488/abundance.tsv",
    )

    normal_samples = concatenate_files(
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/New/ERR358485/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/New/ERR358486/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/Old/ERR358485/abundance.tsv",
        f"{homepath}/RNA-seq/RNA-seq-breast-cancer/Kallisto/Output/Old/ERR358486/abundance.tsv",
    )