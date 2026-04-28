"""
02_gc_content_sliding_window.py

Educational bioinformatics demo for:
- Reading DNA FASTA files
- Calculating GC-content profiles
- Applying triangular weighted sliding-window smoothing
- Plotting local GC-content variation
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np


def read_single_fasta(file_path: str | Path) -> Tuple[str, str]:
    """
    Read the first sequence from a FASTA file.

    Args:
        file_path: Path to FASTA file.

    Returns:
        Tuple of sequence and sequence name.
    """
    sequence_parts = []
    name = "unnamed_sequence"

    with open(file_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                name = line[1:].strip() or "unnamed_sequence"
            else:
                sequence_parts.append(line.upper().replace("U", "T"))

    sequence = clean_dna_sequence("".join(sequence_parts))

    if not sequence:
        raise ValueError(f"No valid DNA sequence found in {file_path}")

    return sequence, name


def clean_dna_sequence(sequence: str) -> str:
    """
    Keep only A, T, G, and C.

    Args:
        sequence: Raw DNA sequence.

    Returns:
        Cleaned DNA sequence.
    """
    valid_bases = {"A", "T", "G", "C"}
    return "".join(base for base in sequence.upper().replace("U", "T") if base in valid_bases)


def triangular_weights(window_size: int) -> np.ndarray:
    """
    Create triangular weights for weighted moving average.

    Example for window_size=9:
    [1, 2, 3, 4, 5, 4, 3, 2, 1]

    Args:
        window_size: Odd window size.

    Returns:
        Array of triangular weights.
    """
    if window_size % 2 == 0:
        raise ValueError("window_size must be odd for triangular weighting.")

    half_window = window_size // 2
    increasing = list(range(1, half_window + 2))
    decreasing = list(range(half_window, 0, -1))

    return np.array(increasing + decreasing)


def gc_content_index(sequence: str, window_size: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate weighted GC-content profile.

    Args:
        sequence: DNA sequence.
        window_size: Odd sliding-window size.

    Returns:
        Positions and weighted GC-content values.
    """
    sequence = clean_dna_sequence(sequence)

    if window_size > len(sequence):
        raise ValueError("window_size cannot be larger than sequence length.")

    weights = triangular_weights(window_size)
    half_window = window_size // 2

    base_values = np.array([
        1.0 if base in {"G", "C"} else 0.0
        for base in sequence
    ])

    x_data = []
    y_data = []

    for center in range(half_window, len(sequence) - half_window):
        start = center - half_window
        end = center + half_window + 1

        window_values = base_values[start:end]
        weighted_value = np.sum(window_values * weights) / np.sum(weights)

        x_data.append(center + 1)
        y_data.append(weighted_value)

    return np.array(x_data), np.array(y_data)


def plot_gc_content(
    x_data: np.ndarray,
    y_data: np.ndarray,
    sequence_length: int,
    window_size: int,
    sequence_name: str,
    output_file: str | Path,
) -> None:
    """
    Plot GC-content profile.

    Args:
        x_data: Sequence positions.
        y_data: GC-content values.
        sequence_length: Total sequence length.
        window_size: Window size.
        sequence_name: Name of sequence.
        output_file: Output PNG path.
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 4))
    plt.plot(x_data, y_data, linewidth=1.5)
    plt.axhline(y=0.5, linestyle="--", linewidth=1)
    plt.xlim(1, sequence_length)
    plt.ylim(0, 1)
    plt.xlabel("Base number")
    plt.ylabel(f"GC content weighted moving average over {window_size} bases")
    plt.title(sequence_name)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def run_demo(window_size: int = 21) -> Tuple[np.ndarray, np.ndarray, int, str]:
    """
    Run built-in demo sequence.

    Args:
        window_size: Odd window size.

    Returns:
        x data, y data, sequence length, sequence name.
    """
    sequence_name = "demo_dna_sequence"
    sequence = (
        "ATATATATATATATATATAT"
        "GCGCGCGCGCGCGCGCGCGC"
        "ATATATATATATATATATAT"
        "GGGGCCCCGGGGCCCCGGGG"
    )

    x_data, y_data = gc_content_index(sequence, window_size)

    return x_data, y_data, len(sequence), sequence_name


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate and plot weighted GC-content profile."
    )

    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Input DNA FASTA file. If omitted, a demo sequence is used.",
    )

    parser.add_argument(
        "--window-size",
        type=int,
        default=101,
        help="Odd sliding-window size.",
    )

    parser.add_argument(
        "--output",
        type=str,
        default="figures/gc_content_profile.png",
        help="Output plot path.",
    )

    args = parser.parse_args()

    if args.input is None:
        print("No input file provided. Running built-in demo.")
        x_data, y_data, sequence_length, sequence_name = run_demo(window_size=21)
        window_size = 21
    else:
        sequence, sequence_name = read_single_fasta(args.input)
        x_data, y_data = gc_content_index(sequence, args.window_size)
        sequence_length = len(sequence)
        window_size = args.window_size

    plot_gc_content(
        x_data=x_data,
        y_data=y_data,
        sequence_length=sequence_length,
        window_size=window_size,
        sequence_name=sequence_name,
        output_file=args.output,
    )

    print(f"GC-content plot saved to: {args.output}")


if __name__ == "__main__":
    main()
