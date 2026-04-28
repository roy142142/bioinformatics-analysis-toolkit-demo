"""
Educational bioinformatics demo for:
- Reading multiple DNA sequences from a FASTA file
- Trimming sequences to equal length
- Splitting sequences into fixed-size windows
- Plotting sequence-logo style visualizations using Biotite
Note:
- This script requires the `biotite` package.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.graphics as graphics


def read_fasta_records(file_path: str | Path) -> Dict[str, str]:
    """
    Read DNA FASTA records.
    """
    records: Dict[str, List[str]] = {}
    current_name: str | None = None

    with open(file_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                current_name = line[1:].strip() or "unnamed_record"
                records[current_name] = []
            else:
                if current_name is None:
                    raise ValueError("FASTA file must start with a FASTA header.")
                records[current_name].append(line.upper().replace("U", "T"))

    return {
        name: clean_dna_sequence("".join(parts))
        for name, parts in records.items()
    }


def clean_dna_sequence(sequence: str) -> str:
    """
    Keep only A, T, G, and C.
    """
    valid_bases = {"A", "T", "G", "C"}
    return "".join(base for base in sequence.upper().replace("U", "T") if base in valid_bases)


def trim_sequences_to_equal_length(records: Dict[str, str]) -> Dict[str, str]:
    """
    Trim all sequences to the length of the shortest sequence.

    Args:
        records: Dictionary of DNA sequences.

    Returns:
        Trimmed sequence dictionary.
    """
    if not records:
        raise ValueError("No sequences provided.")

    min_length = min(len(sequence) for sequence in records.values())

    if min_length == 0:
        raise ValueError("At least one sequence is empty.")

    return {
        name: sequence[:min_length]
        for name, sequence in records.items()
    }


def get_sequence_window(
    records: Dict[str, str],
    window_index: int = 0,
    window_size: int = 100,
    max_sequences: int | None = None,
) -> List[seq.NucleotideSequence]:
    """
    Extract one fixed-size window from all sequences.
    """
    trimmed_records = trim_sequences_to_equal_length(records)

    start = window_index * window_size
    end = start + window_size

    nucleotide_sequences: List[seq.NucleotideSequence] = []

    for index, sequence_string in enumerate(trimmed_records.values()):
        if max_sequences is not None and index >= max_sequences:
            break

        if end > len(sequence_string):
            raise ValueError(
                f"Window {window_index} with size {window_size} exceeds sequence length."
            )

        nucleotide_sequences.append(
            seq.NucleotideSequence(sequence_string[start:end])
        )

    return nucleotide_sequences


def create_identity_alignment(
    sequences: List[seq.NucleotideSequence],
) -> align.Alignment:
    """
    Create an alignment object assuming sequences are already equal length
    and positionally aligned.
    """
    if not sequences:
        raise ValueError("No sequences provided for alignment.")

    sequence_length = len(sequences[0])

    if any(len(sequence) != sequence_length for sequence in sequences):
        raise ValueError("All sequences must have the same length.")

    trace = (
        np.tile(np.arange(sequence_length), len(sequences))
        .reshape(len(sequences), sequence_length)
        .transpose()
    )

    return align.Alignment(sequences=sequences, trace=trace)


def plot_sequence_logo(
    alignment: align.Alignment,
    output_file: str | Path,
    figure_width: float = 20.0,
    figure_height: float = 5.0,
) -> None:
    """
    Plot sequence logo.

    Args:
        alignment: Biotite alignment.
        output_file: Output PNG file.
        figure_width: Figure width.
        figure_height: Figure height.
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(figure_width, figure_height))
    ax = fig.add_subplot(111)

    graphics.plot_sequence_logo(ax, alignment)
    ax.axis("off")

    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate sequence-logo visualization from DNA FASTA records."
    )

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input multi-record DNA FASTA file.",
    )

    parser.add_argument(
        "--window-index",
        type=int,
        default=0,
        help="Zero-based window index to visualize.",
    )

    parser.add_argument(
        "--window-size",
        type=int,
        default=100,
        help="Window size in bases.",
    )

    parser.add_argument(
        "--max-sequences",
        type=int,
        default=None,
        help="Optional maximum number of sequences to include.",
    )

    parser.add_argument(
        "--output",
        type=str,
        default="figures/sequence_logo.png",
        help="Output sequence-logo PNG file.",
    )

    args = parser.parse_args()

    records = read_fasta_records(args.input)

    sequences = get_sequence_window(
        records=records,
        window_index=args.window_index,
        window_size=args.window_size,
        max_sequences=args.max_sequences,
    )

    alignment = create_identity_alignment(sequences)

    plot_sequence_logo(
        alignment=alignment,
        output_file=args.output,
    )

    print(f"Sequence logo saved to: {args.output}")


if __name__ == "__main__":
    main()
