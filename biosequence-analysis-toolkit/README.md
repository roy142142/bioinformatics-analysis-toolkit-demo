
# Biosequence Analysis Toolkit

## Overview

This module contains an educational Python-based toolkit for DNA and protein sequence analysis. The project was originally developed as part of an advanced programming course in biotechnology and has been cleaned, reorganized, and documented for portfolio use.

The toolkit demonstrates how biological sequences can be processed computationally to extract useful sequence-level and protein-level features.

## Project Motivation

DNA sequences encode protein sequences, and protein sequences determine many structural and functional properties of biological molecules. This module demonstrates a basic workflow for moving from nucleotide sequence analysis to protein feature estimation.

The module includes two general analysis directions:

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

## Module Structure

```text
biosequence-analysis-toolkit/
│
├── README.md
├── scripts/
│   ├── 01_orf_translation.py
│   ├── 02_gc_content_sliding_window.py
│   ├── 03_protein_properties.py
│   └── 04_hydropathy_plot.py
│
├── data/
│   └── README.md
│
└── figures/
    └── .gitkeep
