# Biomedical Web Scraping Demo

## Overview

This module contains an educational Python web scraping example for collecting structured biomedical information from publicly available web pages.

The original version of this project was developed as a practical exercise in extracting biomedical and pharmaceutical information from HTML pages. The code has been cleaned and reorganized for portfolio use, with the aim of demonstrating general web scraping skills rather than providing a production database.

The module demonstrates how to:

- Generate multiple index-page URLs programmatically
- Send HTTP requests using Python
- Parse HTML content
- Extract text and links using XPath selectors
- Collect structured records across multiple pages
- Save extracted information into a CSV file
- Handle missing pages and failed requests safely

## Project Motivation

Biomedical information is often distributed across public web resources. Programmatic data collection can be useful for building small reference datasets, preparing inputs for downstream analysis, and practicing data extraction workflows.

This project demonstrates the basic logic behind biomedical web scraping and structured data collection.

## Module Structure

```text
biomedical-web-scraping-demo/
│
├── README.md
├── scripts/
    └── biomedical_index_scraper_demo.py
