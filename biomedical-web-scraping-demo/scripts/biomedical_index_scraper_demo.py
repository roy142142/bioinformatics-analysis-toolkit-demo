"""
biomedical_index_scraper_demo.py

Educational biomedical web scraping demo.

This script demonstrates how to:
- Generate alphabetical index-page URLs
- Send HTTP requests
- Parse HTML pages
- Extract item names and links using XPath
- Save structured results into a CSV file

Important:
This script is intended for educational and portfolio purposes only.
"""

from __future__ import annotations

import csv
import string
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import requests
from lxml import html

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

@dataclass
class ScraperConfig:
    """
    Configuration for the scraping workflow.
    """

    base_url: str = "https://example.org"
    index_url_pattern: str = "https://example.org/index/{suffix}.html"

    # Example XPath selector.
    item_xpath: str = "//a/text()"
    link_xpath: str = "//a/@href"

    output_file: str = "biomedical-web-scraping-demo/outputs/example_output.csv"
    request_delay_seconds: float = 1.0
    timeout_seconds: int = 15

    user_agent: str = (
        "Mozilla/5.0 (compatible; EducationalBiomedicalScraper/1.0)"
    )


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def generate_alphabetical_suffixes(include_numbers: bool = True) -> List[str]:
    """
    Generate two-character alphabetical suffixes.
    """
    first_characters = list(string.ascii_lowercase)
    second_characters = list(string.ascii_lowercase)

    if include_numbers:
        second_characters.append("0-9")

    suffixes = []

    for first_character in first_characters:
        for second_character in second_characters:
            suffixes.append(first_character + second_character)

    return suffixes


def fetch_html(url: str, config: ScraperConfig) -> Optional[str]:
    """
    Request a web page and return HTML text.
    """
    headers = {
        "User-Agent": config.user_agent
    }

    try:
        response = requests.get(
            url,
            headers=headers,
            allow_redirects=False,
            timeout=config.timeout_seconds,
        )
    except requests.RequestException as error:
        print(f"Request failed for {url}: {error}")
        return None

    if response.status_code != 200:
        print(f"Skipping {url} | Status code: {response.status_code}")
        return None

    return response.text


def extract_items_from_html(
    html_text: str,
    source_url: str,
    config: ScraperConfig,
) -> List[dict]:
    """
    Extract item names and links from an HTML page.
    """
    tree = html.fromstring(html_text)

    item_names = tree.xpath(config.item_xpath)
    item_links = tree.xpath(config.link_xpath)

    records = []

    # Keep a conservative pairing approach.
    # If the number of extracted names and links differs,
    # this will only pair available matching indices.
    max_records = min(len(item_names), len(item_links))

    for index in range(max_records):
        item_name = str(item_names[index]).strip()
        relative_link = str(item_links[index]).strip()

        if not item_name:
            continue

        records.append(
            {
                "item_name": item_name,
                "relative_link": relative_link,
                "source_url": source_url,
            }
        )

    return records


def save_records_to_csv(records: List[dict], output_file: str | Path) -> None:
    """
    Save records to CSV.
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ["item_name", "relative_link", "source_url"]

    with output_file.open("w", encoding="utf-8", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


# ---------------------------------------------------------------------
# Main scraping workflow
# ---------------------------------------------------------------------

def run_index_scraper(config: ScraperConfig) -> List[dict]:
    """
    Run scraping over generated alphabetical index pages.
    """
    suffixes = generate_alphabetical_suffixes(include_numbers=True)
    all_records = []

    print(f"Total index pages to check: {len(suffixes)}")

    for suffix in suffixes:
        url = config.index_url_pattern.format(suffix=suffix)
        print(f"Checking: {url}")

        html_text = fetch_html(url, config)

        if html_text is None:
            time.sleep(config.request_delay_seconds)
            continue

        page_records = extract_items_from_html(
            html_text=html_text,
            source_url=url,
            config=config,
        )

        print(f"Extracted records: {len(page_records)}")

        all_records.extend(page_records)

        time.sleep(config.request_delay_seconds)

    return all_records


def main() -> None:
    """
    Entry point for the scraping demo.
    """
    config = ScraperConfig()

    records = run_index_scraper(config)

    save_records_to_csv(
        records=records,
        output_file=config.output_file,
    )

    print(f"Total extracted records: {len(records)}")
    print(f"Output saved to: {config.output_file}")


if __name__ == "__main__":
    main()
