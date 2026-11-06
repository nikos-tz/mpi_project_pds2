import os
import sys
import glob
import re
import argparse

import numpy as np
import matplotlib.pyplot as plt


def read_binary_file(filename):
    """
    Read a binary dataset file.

    Format:
        int64: D
        int64: N
        double[N][D]: data
    """

    with open(filename, "rb") as fp:

        header = np.fromfile(fp, dtype=np.int64, count=2)

        if len(header) != 2:
            raise ValueError(f"Invalid header in {filename}")

        d, n = header

        data = np.fromfile(fp, dtype=np.float64)

        expected_size = n * d

        if len(data) != expected_size:
            raise ValueError(
                f"Invalid data size in {filename}. "
                f"Expected {expected_size} values, got {len(data)}."
            )

        data = data.reshape((n, d))

    return n, d, data


def main():

    parser = argparse.ArgumentParser(
        description="Plot MPI processes' local data and pivot."
    )

    parser.add_argument(
        "directory",
        help="Directory containing proc_data_rank_*.bin and pivot.bin"
    )

    args = parser.parse_args()

    directory = args.directory

    if not os.path.isdir(directory):
        sys.exit(f"Error: Directory '{directory}' does not exist.")

    # Find process data files
    files = glob.glob(
        os.path.join(directory, "proc_data_rank_*.bin")
    )

    if not files:
        sys.exit("Error: No proc_data_rank_*.bin files found.")

    # Sort files according to rank number
    files.sort(
        key=lambda f: int(
            re.search(r"proc_data_rank_(\d+)\.bin$", f).group(1)
        )
    )

    # Read pivot
    pivot_file = os.path.join(directory, "pivot.bin")

    if not os.path.isfile(pivot_file):
        sys.exit("Error: pivot.bin not found.")

    pivot_n, pivot_d, pivot = read_binary_file(pivot_file)

    if pivot_n != 1:
        sys.exit("Error: pivot.bin must contain exactly one point.")

    # Check dimensionality
    if pivot_d != 2:
        sys.exit(
            f"Cannot plot the dataset: d = {pivot_d}. "
            "This script supports only 2-dimensional points."
        )

    # Read data from each process
    process_data = {}

    for filename in files:

        match = re.search(
            r"proc_data_rank_(\d+)\.bin$", filename
        )

        rank = int(match.group(1))

        n, d, data = read_binary_file(filename)

        if d != pivot_d:
            sys.exit(
                f"Error: Rank {rank} has d = {d}, "
                f"but pivot has d = {pivot_d}."
            )

        process_data[rank] = data

        print(
            f"Rank {rank}: {n} points loaded."
        )

    # Plot
    plt.figure(figsize=(10, 8))

    for rank, data in process_data.items():

        plt.scatter(
            data[:, 0],
            data[:, 1],
            label=f"Rank {rank}",
            alpha=0.6,
            s=15
        )

    # Plot pivot
    plt.scatter(
        pivot[0, 0],
        pivot[0, 1],
        color="red",
        marker="x",
        s=150,
        linewidths=3,
        label="Pivot",
        zorder=10
    )

    plt.title("Distribution of Points Among MPI Processes")

    plt.xlabel("Dimension 1")
    plt.ylabel("Dimension 2")

    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()