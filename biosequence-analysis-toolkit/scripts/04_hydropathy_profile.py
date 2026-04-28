"""
04_hydropathy_profile.py

Educational bioinformatics demo for:
- Reading protein FASTA files
- Calculating Kyte-Doolittle hydropathy values
- Applying triangular weighted moving average
- Plotting hydropathy profiles
- Optionally highlighting manually provided protein regions

"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np


KYTE_DOOLITTLE_SCALE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def read_single_fasta(file_path: str | Path) -> Tuple[str, str]:
    """
    Read the first protein sequence from a FASTA file.
    """
    sequence_parts = []
    name = "unnamed_protein"

    with open(file_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                name = line[1:].strip() or "unnamed_protein"
            else:
                sequence_parts.append(line.upper().replace("*", ""))

    sequence = clean_protein_sequence("".join(sequence_parts))

    if not sequence:
        raise ValueError(f"No valid protein sequence found in {file_path}")

    return sequence, name


def clean_protein_sequence(sequence: str) -> str:
    """
    Keep only amino acids supported by the Kyte-Doolittle scale.
    """
    return "".join(
        aa for aa in sequence.upper().replace("*", "")
        if aa in KYTE_DOOLITTLE_SCALE
    )


def triangular_weights(window_size: int) -> np.ndarray:
    """
    Create triangular weights.
    """
    if window_size % 2 == 0:
        raise ValueError("window_size must be odd for triangular weighting.")

    half_window = window_size // 2
    increasing = list(range(1, half_window + 2))
    decreasing = list(range(half_window, 0, -1))

    return np.array(increasing + decreasing)


def hydropathy_index(sequence: str, window_size: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate weighted hydropathy profile.
    """
    sequence = clean_protein_sequence(sequence)

    if window_size > len(sequence):
        raise ValueError("window_size cannot be larger than protein length.")

    values = np.array([
        KYTE_DOOLITTLE_SCALE[aa]
        for aa in sequence
    ])

    weights = triangular_weights(window_size)
    half_window = window_size // 2

    x_data = []
    y_data = []

    for center in range(half_window, len(sequence) - half_window):
        start = center - half_window
        end = center + half_window + 1

        window_values = values[start:end]
        weighted_value = np.sum(window_values * weights) / np.sum(weights)

        x_data.append(center + 1)
        y_data.append(weighted_value)

    return np.array(x_data), np.array(y_data)


def identify_hydrophobic_regions(
    x_data: np.ndarray,
    y_data: np.ndarray,
    threshold: float = 1.6,
) -> List[Tuple[int, int]]:
    """
    Identify continuous regions above the hydrophobicity threshold.
    """
    regions: List[Tuple[int, int]] = []
    in_region = False
    region_start = None

    for position, value in zip(x_data, y_data):
        if value >= threshold and not in_region:
            in_region = True
            region_start = int(position)

        elif value < threshold and in_region:
            in_region = False
            regions.append((int(region_start), int(position - 1)))
            region_start = None

    if in_region and region_start is not None:
        regions.append((int(region_start), int(x_data[-1])))

    return regions


def plot_hydropathy_profile(
    x_data: np.ndarray,
    y_data: np.ndarray,
    protein_length: int,
    window_size: int,
    protein_name: str,
    output_file: str | Path,
    threshold: float = 1.6,
    highlight_regions: List[Tuple[int, int, str]] | None = None,
) -> None:

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 4))

    if highlight_regions:
        for start, end, region_type in highlight_regions:
            if region_type.lower() == "helix":
                color = "pink"
            elif region_type.lower() == "sheet":
                color = "orange"
            else:
                color = "lightgrey"

            plt.axvspan(start, end, facecolor=color, alpha=0.5)

    plt.plot(x_data, y_data, linewidth=1.8)
    plt.axhline(y=threshold, linestyle="--", linewidth=1)
    plt.axhline(y=0, linestyle=":", linewidth=1)
    plt.xlim(1, protein_length)
    plt.xlabel("Residue number")
    plt.ylabel(f"Hydropathy weighted moving average over {window_size} residues")
    plt.title(protein_name)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def run_demo(window_size: int = 9) -> Tuple[np.ndarray, np.ndarray, int, str]:
    protein_name = "demo_membrane_like_protein"
    protein_sequence = (
        "MKKLLPTAAAGLLLLAAQPAMA"
        "VVVVVVVVLLLLLLLLIIIIII"
        "GSGSGSGSGS"
        "DDDDKKKKRRRR"
    )

    x_data, y_data = hydropathy_index(protein_sequence, window_size)

    return x_data, y_data, len(protein_sequence), protein_name


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate and plot protein hydropathy profile."
    )

    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Input protein FASTA file. If omitted, a demo sequence is used.",
    )

    parser.add_argument(
        "--window-size",
        type=int,
        default=19,
        help="Odd sliding-window size.",
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=1.6,
        help="Hydrophobicity threshold.",
    )

    parser.add_argument(
        "--output",
        type=str,
        default="figures/hydropathy_profile.png",
        help="Output plot path.",
    )

    args = parser.parse_args()

    if args.input is None:
        print("No input file provided. Running built-in demo.")
        x_data, y_data, protein_length, protein_name = run_demo(window_size=9)
        window_size = 9
    else:
        protein_sequence, protein_name = read_single_fasta(args.input)
        x_data, y_data = hydropathy_index(protein_sequence, args.window_size)
        protein_length = len(protein_sequence)
        window_size = args.window_size

    plot_hydropathy_profile(
        x_data=x_data,
        y_data=y_data,
        protein_length=protein_length,
        window_size=window_size,
        protein_name=protein_name,
        output_file=args.output,
        threshold=args.threshold,
    )

    regions = identify_hydrophobic_regions(x_data, y_data, threshold=args.threshold)

    print(f"Hydropathy plot saved to: {args.output}")
    print(f"Predicted hydrophobic regions above threshold {args.threshold}: {regions}")


if __name__ == "__main__":
    main()
