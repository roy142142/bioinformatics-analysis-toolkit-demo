"""
03_protein_properties_electrophoresis.py

Educational bioinformatics demo for:
- Translating ORFs from DNA FASTA records
- Estimating protein molecular weight
- Estimating isoelectric point
- Estimating net charge at pH 7
- Generating simple electrophoresis-style scatter plots
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from Bio.SeqUtils.ProtParam import ProteinAnalysis


START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}


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


def find_orfs(sequence: str, frame: int) -> List[str]:
    """
    Find ORFs in one forward reading frame.

    Args:
        sequence: DNA sequence.
        frame: Reading frame, 0/1/2.

    Returns:
        List of ORF DNA sequences.
    """
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


def translate_sequence(dna_sequence: str) -> str:
    """
    Translate DNA sequence into protein sequence.
    """
    return str(Seq(dna_sequence).translate(to_stop=True))


def extract_proteins_from_fasta(
    file_path: str | Path,
    min_protein_length: int = 30,
) -> List[Dict[str, object]]:
    """
    Extract translated proteins from all ORFs in all forward frames.
    """
    records = read_fasta_records(file_path)
    proteins: List[Dict[str, object]] = []

    for record_name, sequence in records.items():
        for frame in range(3):
            orfs = find_orfs(sequence, frame)

            for index, orf in enumerate(orfs, start=1):
                protein_sequence = translate_sequence(orf)

                if len(protein_sequence) < min_protein_length:
                    continue

                proteins.append(
                    {
                        "record_name": record_name,
                        "frame": frame,
                        "orf_index": index,
                        "protein_sequence": protein_sequence,
                    }
                )

    return proteins


def calculate_manual_molecular_weight(protein_sequence: str) -> float:
    """
    Estimate protein molecular weight using residue masses and water loss.
    """
    residue_masses = {
        "A": 89, "V": 117, "L": 131, "I": 131, "P": 115,
        "F": 165, "W": 204, "M": 149, "G": 75, "S": 105,
        "C": 121, "T": 119, "Y": 181, "N": 132, "Q": 146,
        "D": 133, "E": 147, "K": 146, "R": 174, "H": 155,
    }

    sequence = protein_sequence.upper().replace("*", "")
    total_mass = sum(residue_masses.get(aa, 0) for aa in sequence)

    if len(sequence) <= 1:
        return float(total_mass)

    water_loss = 18 * (len(sequence) - 1)

    return float(total_mass - water_loss)


def calculate_protein_properties(
    protein_record: Dict[str, object],
    ph: float = 7.0,
) -> Dict[str, object]:
    """
    Calculate molecular weight, pI and charge.
    """
    sequence = str(protein_record["protein_sequence"]).replace("*", "")

    analysis = ProteinAnalysis(sequence)
    isoelectric_point = analysis.isoelectric_point()
    charge_at_ph = IsoelectricPoint(sequence).charge_at_pH(ph)

    return {
        "record_name": protein_record["record_name"],
        "frame": protein_record["frame"],
        "orf_index": protein_record["orf_index"],
        "protein_length": len(sequence),
        "manual_molecular_weight_da": round(calculate_manual_molecular_weight(sequence), 3),
        "biopython_molecular_weight_da": round(analysis.molecular_weight(), 3),
        "isoelectric_point": round(isoelectric_point, 3),
        f"charge_at_pH_{ph}": round(charge_at_ph, 3),
        "protein_sequence": sequence,
    }


def write_properties_table(results: List[Dict[str, object]], output_file: str | Path) -> None:
    """
    Save protein properties to CSV.
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if not results:
        raise ValueError("No protein property results to write.")

    with output_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)


def plot_electrophoresis_style(
    results: List[Dict[str, object]],
    output_file: str | Path,
) -> None:
    """
    Plot molecular weight against isoelectric point.
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    molecular_weights = [record["manual_molecular_weight_da"] for record in results]
    isoelectric_points = [record["isoelectric_point"] for record in results]

    plt.figure(figsize=(6, 5))
    plt.scatter(molecular_weights, isoelectric_points, marker="s")
    plt.xlabel("Protein molecular weight (Da)")
    plt.ylabel("Isoelectric point")
    plt.title("Electrophoresis-style protein feature plot")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_charge_by_isoelectric_point(
    results: List[Dict[str, object]],
    ph: float,
    output_file: str | Path,
) -> None:
    """
    Plot charge at selected pH against isoelectric point.
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    charges = [record[f"charge_at_pH_{ph}"] for record in results]
    isoelectric_points = [record["isoelectric_point"] for record in results]

    plt.figure(figsize=(6, 5))
    plt.scatter(isoelectric_points, charges, marker="o")
    plt.xlabel("Isoelectric point")
    plt.ylabel(f"Protein charge at pH {ph}")
    plt.title("Protein charge by isoelectric point")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def run_demo(ph: float = 7.0) -> List[Dict[str, object]]:
    """
    Run demo using a few protein sequences.
    """
    demo_proteins = [
        {
            "record_name": "demo_protein_A",
            "frame": 0,
            "orf_index": 1,
            "protein_sequence": "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE",
        },
        {
            "record_name": "demo_protein_B",
            "frame": 1,
            "orf_index": 1,
            "protein_sequence": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGER",
        },
    ]

    return [
        calculate_protein_properties(record, ph=ph)
        for record in demo_proteins
    ]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Estimate protein properties and generate feature plots."
    )

    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Input DNA FASTA file. If omitted, demo proteins are used.",
    )

    parser.add_argument(
        "--min-protein-length",
        type=int,
        default=30,
        help="Minimum translated protein length to keep.",
    )

    parser.add_argument(
        "--ph",
        type=float,
        default=7.0,
        help="pH value for charge estimation.",
    )

    parser.add_argument(
        "--output-table",
        type=str,
        default="outputs/protein_properties.csv",
        help="Output CSV table.",
    )

    parser.add_argument(
        "--output-dir",
        type=str,
        default="figures",
        help="Output directory for plots.",
    )

    args = parser.parse_args()

    if args.input is None:
        print("No input file provided. Running built-in demo proteins.")
        results = run_demo(ph=args.ph)
    else:
        proteins = extract_proteins_from_fasta(
            args.input,
            min_protein_length=args.min_protein_length,
        )
        results = [
            calculate_protein_properties(protein, ph=args.ph)
            for protein in proteins
        ]

    write_properties_table(results, args.output_table)

    output_dir = Path(args.output_dir)

    plot_electrophoresis_style(
        results,
        output_dir / "electrophoresis_style_plot.png",
    )

    plot_charge_by_isoelectric_point(
        results,
        ph=args.ph,
        output_file=output_dir / "charge_by_isoelectric_point.png",
    )

    print(f"Protein property table saved to: {args.output_table}")
    print(f"Plots saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
