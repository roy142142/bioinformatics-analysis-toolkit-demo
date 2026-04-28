"""
01_orf_translation.py

Educational bioinformatics demo for:
- Reading multi-record FASTA files
- Extracting DNA sequences
- Detecting ORFs in three forward reading frames
- Translating ORFs into protein sequences
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

from Bio.Seq import Seq


START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}


def read_fasta_records(file_path: str | Path) -> Dict[str, str]:
    """
    Read sequences from a FASTA file.

    Args:
        file_path: Path to FASTA file.

    Returns:
        Dictionary mapping record names to DNA sequences.
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
                    raise ValueError("FASTA file must start with a header line beginning with '>'.")
                records[current_name].append(line.upper().replace("U", "T"))

    return {
        name: clean_dna_sequence("".join(sequence_lines))
        for name, sequence_lines in records.items()
    }


def clean_dna_sequence(sequence: str) -> str:
    """
    Keep only standard DNA bases.

    Args:
        sequence: Raw DNA sequence.

    Returns:
        Cleaned DNA sequence.
    """
    valid_bases = {"A", "T", "G", "C"}
    return "".join(base for base in sequence.upper().replace("U", "T") if base in valid_bases)


def find_orfs(sequence: str, frame: int) -> List[str]:
    """
    Find ORFs in one forward reading frame.

    Args:
        sequence: DNA sequence.
        frame: Reading frame. Must be 0, 1, or 2.

    Returns:
        List of ORF DNA sequences, including stop codons.
    """
    if frame not in {0, 1, 2}:
        raise ValueError("frame must be 0, 1, or 2.")

    sequence = clean_dna_sequence(sequence)
    orfs: List[str] = []

    for start in range(frame, len(sequence) - 2, 3):
        codon = sequence[start:start + 3]

        if codon != START_CODON:
            continue

        for stop in range(start + 3, len(sequence) - 2, 3):
            stop_codon = sequence[stop:stop + 3]

            if stop_codon in STOP_CODONS:
                orfs.append(sequence[start:stop + 3])
                break

    return orfs


def find_orfs_all_frames(sequence: str) -> Dict[int, List[str]]:
    """
    Find ORFs in three forward reading frames.

    Args:
        sequence: DNA sequence.

    Returns:
        Dictionary mapping frame number to ORF list.
    """
    return {
        frame: find_orfs(sequence, frame)
        for frame in range(3)
    }


def translate_orfs(orfs: List[str]) -> List[str]:
    """
    Translate DNA ORFs into amino-acid sequences.

    Args:
        orfs: List of DNA ORF sequences.

    Returns:
        List of translated protein sequences.
    """
    proteins: List[str] = []

    for orf in orfs:
        protein = str(Seq(orf).translate(to_stop=True))
        proteins.append(protein)

    return proteins


def analyze_fasta_file(file_path: str | Path, min_protein_length: int = 30) -> List[Dict[str, object]]:
    """
    Run ORF detection and translation for all records in a FASTA file.

    Args:
        file_path: Input FASTA file.
        min_protein_length: Minimum amino-acid length to keep.

    Returns:
        List of ORF/protein records.
    """
    records = read_fasta_records(file_path)
    results: List[Dict[str, object]] = []

    for record_name, sequence in records.items():
        orfs_by_frame = find_orfs_all_frames(sequence)

        for frame, orfs in orfs_by_frame.items():
            proteins = translate_orfs(orfs)

            for index, protein in enumerate(proteins, start=1):
                if len(protein) < min_protein_length:
                    continue

                results.append(
                    {
                        "record_name": record_name,
                        "frame": frame,
                        "orf_index": index,
                        "protein_length": len(protein),
                        "protein_sequence": protein,
                    }
                )

    return results


def write_translation_results(results: List[Dict[str, object]], output_file: str | Path) -> None:
    """
    Write translated ORF results to TSV.

    Args:
        results: ORF translation results.
        output_file: Output TSV path.
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with output_file.open("w", encoding="utf-8") as handle:
        handle.write("record_name\tframe\torf_index\tprotein_length\tprotein_sequence\n")

        for result in results:
            handle.write(
                f"{result['record_name']}\t"
                f"{result['frame']}\t"
                f"{result['orf_index']}\t"
                f"{result['protein_length']}\t"
                f"{result['protein_sequence']}\n"
            )


def run_demo() -> List[Dict[str, object]]:
    """
    Run a small built-in demo without external files.

    Returns:
        Demo ORF translation results.
    """
    demo_sequence = (
        "AAACCCATGGCTGCTGCTAAATTTGGGATGAAACCCGGGTTTAAACCC"
        "ATGCCCTTTGGGAAATAGGGG"
    )

    records = {"demo_sequence": demo_sequence}
    results: List[Dict[str, object]] = []

    for record_name, sequence in records.items():
        orfs_by_frame = find_orfs_all_frames(sequence)

        for frame, orfs in orfs_by_frame.items():
            proteins = translate_orfs(orfs)

            for index, protein in enumerate(proteins, start=1):
                results.append(
                    {
                        "record_name": record_name,
                        "frame": frame,
                        "orf_index": index,
                        "protein_length": len(protein),
                        "protein_sequence": protein,
                    }
                )

    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Detect ORFs and translate DNA sequences."
    )

    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Input DNA FASTA file. If omitted, a demo sequence is used.",
    )

    parser.add_argument(
        "--min-protein-length",
        type=int,
        default=30,
        help="Minimum translated protein length to keep.",
    )

    parser.add_argument(
        "--output",
        type=str,
        default="outputs/orf_translation_results.tsv",
        help="Output TSV file.",
    )

    args = parser.parse_args()

    if args.input is None:
        print("No input file provided. Running built-in demo.")
        results = run_demo()
    else:
        results = analyze_fasta_file(
            file_path=args.input,
            min_protein_length=args.min_protein_length,
        )

    write_translation_results(results, args.output)

    print(f"Number of translated ORFs: {len(results)}")
    print(f"Results saved to: {args.output}")


if __name__ == "__main__":
    main()
