# Bioinformatics Analysis Toolkit Demo

Educational Python toolkit for bioinformatics and biomedical data analysis, including DNA/protein sequence analysis, ORF detection, DNA-to-protein translation, GC-content profiling, protein feature estimation, hydropathy analysis, sequence-logo visualization, and biomedical web scraping.

## Overview

This repository contains educational and portfolio-oriented bioinformatics and biomedical data analysis tools developed in Python. The repository is organized as a collection of small analysis modules, each focusing on a specific area of computational biology or biomedical data processing.

The repository currently includes two modules:

1. **Biosequence Analysis Toolkit**  
   A DNA and protein sequence analysis module covering ORF detection, DNA-to-protein translation, GC-content profiling, protein physicochemical feature estimation, hydropathy analysis, and sequence-logo visualization.

2. **Biomedical Web Scraping Demo**  
   A general web scraping module demonstrating how biomedical information can be collected from public web pages, parsed from HTML, and saved into structured files.

## Repository Map

```text
bioinformatics-analysis-toolkit-demo/
│
├── README.md
│
├── biosequence-analysis-toolkit/
│   ├── README.md
│   │
│   ├── scripts/
│   │   ├── 01_orf_translation.py
│   │   ├── 02_gc_content_sliding_window.py
│   │   ├── 03_protein_properties_electrophoresis.py
│   │   ├── 04_hydropathy_profile.py
│   │   └── 05_sequence_logo_alignment.py
│   │
│   └── figures/
│       ├── gc_content_profile_demo.png
│       ├── hydropathy_profile_demo.png
│       └── sequence_logo_demo.png
│
└── biomedical-web-scraping-demo/
    ├── README.md
    │
    ├── scripts/
    └── biomedical_index_scraper_demo.py
