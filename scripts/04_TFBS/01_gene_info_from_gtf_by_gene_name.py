#!/usr/bin/env python3
"""
Step 01: Extract zebrafish gene information from a GTF file.

Inputs:
    data/BET_posi.txt or data/FLU_posi.txt
    data/Danio_rerio.GRCz11.115.gtf

Outputs:
    results/04_TFBS/BET_posi_gene_table.tsv
    results/04_TFBS/FLU_posi_gene_table.tsv

Reference genome: GRCz11
"""

import argparse
import re
from pathlib import Path

import pandas as pd


# --------------------------------------------------
# 1. Define project directories
# --------------------------------------------------

PROJECT_DIR = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "results" / "04_TFBS"

GTF_FILE = DATA_DIR / "Danio_rerio.GRCz11.115.gtf"


# --------------------------------------------------
# 2. Read command-line arguments
# --------------------------------------------------

parser = argparse.ArgumentParser(
    description="Extract zebrafish gene annotations from GTF."
)

parser.add_argument(
    "--compound",
    required=True,
    choices=["BET", "FLU"],
    help="Compound name: BET or FLU"
)

args = parser.parse_args()

compound = args.compound

gene_list_file = DATA_DIR / f"{compound}_posi.txt"
output_file = OUTPUT_DIR / f"{compound}_posi_gene_table.tsv"


# --------------------------------------------------
# 3. Check input files
# --------------------------------------------------

if not GTF_FILE.is_file():
    raise FileNotFoundError(f"GTF file not found: {GTF_FILE}")

if not gene_list_file.is_file():
    raise FileNotFoundError(
        f"Gene list file not found: {gene_list_file}"
    )

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# 4. Load target gene identifiers
# --------------------------------------------------

with open(gene_list_file, encoding="utf-8") as f:
    genes = {
        line.strip()
        for line in f
        if line.strip()
    }

print(f"Compound: {compound}")
print(f"Loaded {len(genes)} unique gene identifiers")


# --------------------------------------------------
# 5. Read GTF gene annotations
# --------------------------------------------------

gtf_columns = [
    "chrom", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute"
]

gtf = pd.read_csv(
    GTF_FILE,
    sep="\t",
    comment="#",
    header=None,
    names=gtf_columns,
    dtype={"chrom": str}
)

# Retain only gene-level annotations.
gtf = gtf.loc[gtf["feature"] == "gene"].copy()


# --------------------------------------------------
# 6. Extract gene attributes
# --------------------------------------------------

def extract_attribute(attribute_string, key):
    """Extract a specified attribute from a GTF record."""
    pattern = rf'(?:^|;\s*){re.escape(key)}\s+"([^"]+)"'
    match = re.search(pattern, str(attribute_string))
    return match.group(1) if match else None


gtf["gene_id"] = gtf["attribute"].apply(
    lambda x: extract_attribute(x, "gene_id")
)

gtf["gene_name"] = gtf["attribute"].apply(
    lambda x: extract_attribute(x, "gene_name")
)

gtf["gene_biotype"] = gtf["attribute"].apply(
    lambda x: (
        extract_attribute(x, "gene_biotype")
        or extract_attribute(x, "gene_type")
        or ""
    )
)


# --------------------------------------------------
# 7. Match target genes to GTF annotations
# --------------------------------------------------

matched = (
    gtf["gene_name"].isin(genes)
    | gtf["gene_id"].isin(genes)
)

selected = gtf.loc[matched].copy()


# --------------------------------------------------
# 8. Prepare output table
# --------------------------------------------------

output_columns = [
    "gene_id",
    "chrom",
    "start",
    "end",
    "strand",
    "gene_name",
    "gene_biotype"
]

output_df = selected[output_columns].copy()

# Standardize chromosome names.
output_df["chrom"] = output_df["chrom"].apply(
    lambda x: x if str(x).startswith("chr") else "chr" + str(x)
)

# Preserve original GTF coordinates in this step.
# GTF uses 1-based inclusive coordinates.
# BED coordinate conversion will be handled in Step 02.

output_df = output_df.drop_duplicates()

output_df.to_csv(
    output_file,
    sep="\t",
    index=False
)


# --------------------------------------------------
# 9. Report matching statistics
# --------------------------------------------------

matched_ids = set(output_df["gene_id"].dropna())
matched_names = set(output_df["gene_name"].dropna())

matched_inputs = genes & (matched_ids | matched_names)
unmatched_inputs = genes - matched_inputs

unmatched_file = OUTPUT_DIR / f"{compound}_unmatched_genes.txt"

with open(unmatched_file, "w", encoding="utf-8") as f:
    for gene in sorted(unmatched_inputs):
        f.write(f"{gene}\n")

print("-" * 50)
print(f"Input gene identifiers: {len(genes)}")
print(f"Matched input identifiers: {len(matched_inputs)}")
print(f"Matched gene records: {len(output_df)}")
print(f"Unmatched identifiers: {len(unmatched_inputs)}")
print(f"Gene table saved to: {output_file}")
print(f"Unmatched gene list saved to: {unmatched_file}")
print("-" * 50)
