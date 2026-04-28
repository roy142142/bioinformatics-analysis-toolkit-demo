# Biosequence Analysis Toolkit

## Overview

This module contains an educational Python-based toolkit for DNA and protein sequence analysis. The project was originally developed as part of an advanced programming course in biotechnology and has been cleaned, reorganized, and documented for portfolio use.

The toolkit demonstrates how biological sequences can be processed computationally to extract useful sequence-level and protein-level features, including ORF detection, DNA-to-protein translation, GC-content profiling, protein molecular weight estimation, isoelectric point estimation, charge estimation, hydropathy profiling, and sequence-logo visualization.

## Project Motivation

DNA sequences encode protein sequences, and protein sequences determine many structural and functional properties of biological molecules. This module demonstrates a basic workflow for moving from nucleotide sequence analysis to protein feature estimation and sequence-level visualization.

The module includes three general analysis directions:

1. **DNA-based analysis**
   - Reading DNA sequences from FASTA-like files
   - Detecting open reading frames in different reading frames
   - Translating DNA sequences into protein sequences
   - Calculating GC-content profiles

2. **Protein-based analysis**
   - Estimating protein molecular weight
   - Estimating isoelectric point and net charge
   - Generating electrophoresis-style feature summaries
   - Calculating and visualizing hydropathy profiles

3. **Comparative sequence visualization**
   - Trimming multiple sequences to equal length
   - Splitting sequences into fixed-size windows
   - Generating sequence-logo style visualizations

## Module Structure

```text
biosequence-analysis-toolkit/
│
├── README.md
├── scripts/
│   ├── 01_orf_translation.py
│   ├── 02_gc_content_sliding_window.py
│   ├── 03_protein_properties_electrophoresis.py
│   ├── 04_hydropathy_profile.py
│   └── 05_sequence_logo_alignment.py
│
│
└── figures/

