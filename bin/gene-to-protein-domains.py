#!/usr/bin/env python3

"""
Retrieve protein domain positions for a gene and output as TSV.

This script only works for human genes (Homo sapiens, taxonomy_id:9606, reviewed entries only).
"""

import sys
import os
import json
import urllib.request
import argparse
import logging
import time
from urllib.error import URLError, HTTPError

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S"
)

def fetch_json_with_retry(url, headers=None, retries=3, delay=2):
    """
    Fetch JSON data from a URL with retry logic.

    Args:
        url (str): The URL to fetch.
        headers (dict, optional): HTTP headers.
        retries (int): Number of retries.
        delay (int): Delay between retries in seconds.

    Returns:
        dict: Parsed JSON response.

    Raises:
        Exception: If all retries fail.
    """
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=headers or {})
            with urllib.request.urlopen(req) as response:
                return json.load(response)
        except (URLError, HTTPError) as e:
            logging.warning(f"Attempt {attempt} failed for {url}: {e}")
            if attempt < retries:
                logging.warning(f"Retrying... (attempt {attempt + 1})")
                time.sleep(delay)
            else:
                raise

class ProteinDomainFetcher:
    """
    Fetch and display protein domain information for a given human gene using UniProt and InterPro APIs.

    NOTE:
        This script only works for human genes (Homo sapiens, taxonomy_id:9606, reviewed entries only).

    Args:
        gene_name (str): The gene name to search for (e.g., 'NOTCH1').
        hide_types (str, optional): Comma-separated string of domain types to hide (e.g., 'family,repeat').

    Example:
        fetcher = ProteinDomainFetcher('NOTCH1', 'family,repeat')
        fetcher.run()
    """
    def __init__(self, gene_name, hide_types=None, retries=3):
        """
        Initialize the fetcher with a gene name and optional types to hide.

        Args:
            gene_name (str): The gene name to search for (e.g., 'NOTCH1').
            hide_types (str, optional): Comma-separated string of domain types to hide (e.g., 'family,repeat').
        """
        self.gene_name = gene_name
        self.hide_types = hide_types
        self.uniprot_id = None
        self.results = None
        self.retries = retries

    def get_uniprot_id(self):
        """
        Retrieve the UniProt ID for the gene name (human, reviewed entries only).

        Returns:
            str or None: UniProt accession string if found, else None.

        Note:
            Only reviewed human entries are considered.

        Example:
            uniprot_id = self.get_uniprot_id()
        """
        url = (
            f"https://rest.uniprot.org/uniprotkb/search?fields=accession&format=json"
            f"&query=((gene:{self.gene_name})+AND+(taxonomy_id:9606)+AND+(reviewed:true))"
        )
        data = fetch_json_with_retry(url, retries=self.retries)
        results = data.get("results", [])
        if not results:
            return None
        return results[0].get("primaryAccession")

    def get_domains(self):
        """
        Retrieve domain information from InterPro for the UniProt ID.

        Returns:
            dict: Parsed JSON response from InterPro API.

        Note:
            Requires self.uniprot_id to be set.

        Example:
            domains = self.get_domains()
        """
        url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/UniProt/{self.uniprot_id}/?page_size=200"
        return fetch_json_with_retry(url, retries=self.retries)

    def print_table(self):
        """
        Print the domain information as a TSV table, optionally hiding specified types. Output is sorted by fragment start position.

        Args:
            None (uses self.results and self.hide_types)

        Note:
            Output is written to a file named <gene_name>_domain.tsv in the current directory.

        Example:
            self.print_table()
        """
        header = [
            "query_gene", "uniprot_id", "pfam_accession", "name", "source_database", "type", "integrated_id", "go_terms",
            "protein_accession", "protein_length", "entry_protein_locations_count",
            "start", "end"
        ]

        hide_types_set = set()
        if self.hide_types:
            hide_types_set = set(t.strip().lower() for t in self.hide_types.split(",") if t.strip())
        rows = set()
        for entry in self.results.get('results', []):
            meta = entry["metadata"]
            typ = meta.get("type", "").lower()
            if typ in hide_types_set:
                continue
            for prot in entry["proteins"]:
                locations = prot.get("entry_protein_locations", [])
                loc_count = len(locations)
                for loc in locations:
                    for frag in loc.get("fragments", []):
                        row = (
                            self.gene_name,
                            self.uniprot_id,
                            meta.get("accession", ""),
                            meta.get("name", ""),
                            meta.get("source_database", ""),
                            meta.get("type", ""),
                            str(meta.get("integrated", "")),
                            str(meta.get("go_terms", "")),
                            prot.get("accession", ""),
                            str(prot.get("protein_length", "")),
                            str(loc_count),
                            str(frag.get("start", "")),
                            str(frag.get("end", ""))
                        )
                        rows.add(row)
        try:
            start_idx = header.index("start")
        except ValueError:
            start_idx = len(header) - 2
        sorted_rows = sorted(
            rows,
            key=lambda r: int(r[start_idx]) if str(r[start_idx]).isdigit() else float('inf')
        )
        output_filename = f"{self.gene_name}_domain.tsv"
        with open(output_filename, "w") as out_f:
            out_f.write("\t".join(header) + "\n")
            for row in sorted_rows:
                out_f.write("\t".join(row) + "\n")

    def run(self):
        logging.info(f"Fetching UniProt ID for: {self.gene_name}")
        self.uniprot_id = self.get_uniprot_id()
        if not self.uniprot_id:
            logging.error(f"No UniProt ID found for: {self.gene_name}")
            sys.exit(1)
        logging.info(f"Found UniProt ID: {self.uniprot_id}")
        logging.info(f"Fetching domain information for: {self.uniprot_id}")
        self.results = self.get_domains()
        logging.info(f"Writing domain table to {self.gene_name}_domain.tsv")
        self.print_table()


class TranscriptInfoFetcher:
    """
    Fetch and save transcript information for a given human gene symbol using Ensembl REST API.

    NOTE:
        This class only works for human genes (Homo sapiens).

    Args:
        gene_name (str): The gene symbol to search for (e.g., 'NOTCH1').

    Example:
        fetcher = TranscriptInfoFetcher('NOTCH1')
        fetcher.run()
    """
    ENSEMBL_SERVER = "https://rest.ensembl.org"
    ENSEMBL_HEADERS = {"Content-Type": "application/json"}

    def __init__(self, gene_name, retries=3):
        """
        Initialize the fetcher with a gene symbol.

        Args:
            gene_name (str): The gene symbol to search for (e.g., 'NOTCH1').
        """
        self.gene_name = gene_name
        self.gene_id = None
        self.gene_info = None
        self.canonical_transcript = None
        self.retries = retries

    def get_gene_id(self):
        """
        Retrieve the Ensembl gene ID for the gene symbol (human).

        Returns:
            str or None: Ensembl gene ID if found, else None.

        Example:
            gene_id = self.get_gene_id()
        """
        url = f"{self.ENSEMBL_SERVER}/lookup/symbol/homo_sapiens/{self.gene_name}"
        try:
            data = fetch_json_with_retry(url, headers=self.ENSEMBL_HEADERS, retries=self.retries)
            return data.get("id")
        except Exception:
            return None

    def get_gene_info(self):
        """
        Retrieve gene information (including MANE and canonical transcript) from Ensembl.

        Returns:
            dict or None: Parsed JSON response from Ensembl API.

        Note:
            Requires self.gene_id to be set.

        Example:
            info = self.get_gene_info()
        """
        url = f"{self.ENSEMBL_SERVER}/lookup/id/{self.gene_id}?mane=1"
        try:
            return fetch_json_with_retry(url, headers=self.ENSEMBL_HEADERS, retries=self.retries)
        except Exception:
            return None

    def print_table(self):
        """
        Write the gene and transcript information as a TSV table to a file.

        Args:
            None (uses self.gene_info)

        Note:
            Output is written to a file named <gene_name>_transcript.tsv in the current directory.

        Example:
            self.print_table()
        """
        if not self.gene_info:
            return
        output_filename = f"{self.gene_name}_transcript.tsv"
        with open(output_filename, "w") as out_f:
            for key in sorted(self.gene_info.keys()):
                out_f.write(f"{key}\t{self.gene_info[key]}\n")

    def run(self):
        logging.info(f"Fetching Ensembl gene ID for: {self.gene_name}")
        self.gene_id = self.get_gene_id()
        if not self.gene_id:
            logging.error(f"No Ensembl gene ID found for: {self.gene_name}")
            sys.exit(1)
        logging.info(f"Found Ensembl gene ID: {self.gene_id}")
        logging.info(f"Fetching transcript information for: {self.gene_id}")
        self.gene_info = self.get_gene_info()
        if not self.gene_info:
            logging.error(f"No gene information found for: {self.gene_name}")
            sys.exit(1)
        self.canonical_transcript = self.gene_info.get("canonical_transcript")
        if self.canonical_transcript:
            self.canonical_transcript = self.canonical_transcript.split(".")[0]
            self.gene_info["canonical_transcript"] = self.canonical_transcript
            logging.info(f"Found MANE transcript: {self.gene_info['canonical_transcript']}")
        logging.info(f"Writing transcript table to {self.gene_name}_transcript.tsv")
        self.print_table()


def main():
    """
    Parse command-line arguments and run the ProteinDomainFetcher and TranscriptInfoFetcher.

    Args:
        None (uses sys.argv)

    Example:
        python {script_name} --gene NOTCH1 --hide family,repeat

    Note:
        This script fetches both protein domain information (from UniProt/InterPro)
        and transcript information (from Ensembl) for a given human gene symbol.
        The outputs are written to two TSV files in the current directory:
            <gene_name>_domain.tsv      - Protein domain annotation table
            <gene_name>_transcript.tsv  - Transcript annotation table
    """
    script_name = os.path.basename(sys.argv[0])
    parser = argparse.ArgumentParser(
        description=(
            f"Retrieve protein domain and transcript information for a human gene and output as TSV files.\n"
            f"\n"
            f"NOTE: This script works with human genes (Homo sapiens, taxonomy_id:9606, reviewed entries only).\n"
            f"\n"
            f"Outputs:\n"
            f"  <gene_name>_domain.tsv      - Protein domain annotation table\n"
            f"  <gene_name>_transcript.tsv  - MANE Transcript annotation table\n"
            f"\n"
            f"Examples:\n"
            f"  python {script_name} --gene NOTCH1\n"
            f"  python {script_name} --gene NOTCH1 --hide family,repeat\n"
            f"  python {script_name} --gene BRCA1 --hide domain\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--gene', required=True,
        help="Gene name to search (e.g. NOTCH1, BRCA1). Example: --gene NOTCH1"
    )
    parser.add_argument(
        '--hide', default=None,
        help="Comma-separated values to hide by 'Type' column in the domain table (e.g. family,repeat). Example: --hide family,repeat"
    )
    parser.add_argument(
        '--fetch', choices=['domain', 'transcript', 'both'], default='both',
        help="Which information to fetch: domain, transcript, or both (default: both)"
    )
    parser.add_argument(
        '--retries', type=int, default=3,
        help="Number of times to retry API requests on failure (default: 3)"
    )
    args = parser.parse_args()

    if args.fetch in ('domain', 'both'):
        domain_fetcher = ProteinDomainFetcher(args.gene.upper(), args.hide, retries=args.retries)
        domain_fetcher.run()
    if args.fetch in ('transcript', 'both'):
        transcript_fetcher = TranscriptInfoFetcher(args.gene.upper(), retries=args.retries)
        transcript_fetcher.run()


if __name__ == "__main__":
    main()
