#!/usr/bin/env python3
"""
Step 02: Associate zebrafish genes with ATAC-seq peaks.

Input:
    results/04_TFBS/{compound}_posi_gene_table.tsv
    data/ATAC_GRCz11_peaks_434221.bed

Output:
    results/04_TFBS/{compound}_posi_CREs.tsv

Method:
    Identify ATAC-seq peaks overlapping gene bodies
    using bedtools intersect (-wa -wb).

Reference genome: GRCz11
"""

import argparse
from pathlib import Path

import pandas as pd
from pybedtools import BedTool
from tqdm import tqdm


# --------------------------------------------------
# 1. Define project directories
# --------------------------------------------------

PROJECT_DIR = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "results" / "04_TFBS"

ATAC_FILE = DATA_DIR / "ATAC_GRCz11_peaks_434221.bed"


# --------------------------------------------------
# 2. Read command-line arguments
# --------------------------------------------------

parser = argparse.ArgumentParser(
    description="Associate gene bodies with ATAC-seq peaks."
)

parser.add_argument(
    "--compound",
    required=True,
    choices=["BET", "FLU"],
    help="Compound name: BET or FLU"
)

args = parser.parse_args()

compound = args.compound

gene_file = OUTPUT_DIR / f"{compound}_posi_gene_table.tsv"
output_file = OUTPUT_DIR / f"{compound}_posi_CREs.tsv"


# --------------------------------------------------
# 3. Check input files
# --------------------------------------------------

if not gene_file.is_file():
    raise FileNotFoundError(
        f"Gene annotation table not found: {gene_file}"
    )

if not ATAC_FILE.is_file():
    raise FileNotFoundError(
        f"ATAC-seq BED file not found: {ATAC_FILE}"
    )

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# 4. Load gene annotation table
# --------------------------------------------------

genes = pd.read_csv(gene_file, sep="\t")

required_columns = [
    "gene_id", "chrom", "start", "end",
    "strand", "gene_name", "gene_biotype"
]

missing_columns = [
    col for col in required_columns
    if col not in genes.columns
]

if missing_columns:
    raise ValueError(
        f"Missing columns in gene table: {missing_columns}"
    )

genes["chrom"] = genes["chrom"].astype(str).apply(
    lambda x: x if x.startswith("chr") else "chr" + x
)

print(f"Compound: {compound}")
print(f"Loaded {len(genes)} gene records")


# --------------------------------------------------
# 5. Convert gene coordinates to BED format
# --------------------------------------------------

# GTF: 1-based inclusive coordinates.
# BED: 0-based, half-open coordinates.
# Therefore, BED start = GTF start - 1.
# BED end = GTF end.

gene_bed_list = []

gene_info_dict = {}

for _, row in genes.iterrows():

    gene_id = str(row["gene_id"])

    start_bed = int(row["start"]) - 1
    end_bed = int(row["end"])

    if start_bed < 0 or end_bed <= start_bed:
        raise ValueError(
            f"Invalid gene coordinates: {gene_id}"
        )

    gene_bed_list.append([
        row["chrom"],
        start_bed,
        end_bed,
        gene_id,
        0,
        row["strand"]
    ])

    gene_info_dict[gene_id] = {
        "gene_name": row["gene_name"],
        "gene_biotype": row["gene_biotype"]
    }

if not gene_bed_list:
    raise ValueError("No genes available for CRE matching.")

gene_bed = BedTool(gene_bed_list)


# --------------------------------------------------
# 6. Load ATAC-seq peak regions
# --------------------------------------------------

# Retain the first three BED columns:
# chromosome, start, and end.

cre_bed = BedTool(str(ATAC_FILE)).cut(range(3))


# --------------------------------------------------
# 7. Identify overlapping gene–peak pairs
# --------------------------------------------------

# -wa: Report the original gene interval.
# -wb: Report the overlapping ATAC-seq interval.
# Any positive interval overlap is retained.

print("Matching ATAC-seq peaks to gene bodies...")

intersect = gene_bed.intersect(
    cre_bed,
    wa=True,
    wb=True
)


# --------------------------------------------------
# 8. Extract gene–CRE associations
# --------------------------------------------------

matches = []

for interval in tqdm(
    intersect,
    desc="Processing gene–peak overlaps"
):

    gene_id = interval.name

    start_gene = int(interval.start)
    end_gene = int(interval.end)

    start_peak = int(interval.fields[7])
    end_peak = int(interval.fields[8])

    gene_info = gene_info_dict[gene_id]

    matches.append({
        "gene_id": gene_id,
        "chrom": interval.chrom,
        "start_bed": start_peak,
        "end_bed": end_peak,
        "start_gene": start_gene,
        "end_gene": end_gene,
        "strand": interval.strand,
        "gene_name": gene_info["gene_name"],
        "gene_biotype": gene_info["gene_biotype"]
    })


# --------------------------------------------------
# 9. Save gene–CRE association table
# --------------------------------------------------

output_columns = [
    "gene_id",
    "chrom",
    "start_bed",
    "end_bed",
    "start_gene",
    "end_gene",
    "strand",
    "gene_name",
    "gene_biotype"
]

matches_df = pd.DataFrame(
    matches,
    columns=output_columns
)

matches_df.to_csv(
    output_file,
    sep="\t",
    index=False
)


# --------------------------------------------------
# 10. Report summary statistics
# --------------------------------------------------

matched_genes = matches_df["gene_id"].nunique()

unique_peaks = matches_df[
    ["chrom", "start_bed", "end_bed"]
].drop_duplicates().shape[0]

print("-" * 50)
print(f"Input gene records: {len(genes)}")
print(f"Genes with overlapping ATAC peaks: {matched_genes}")
print(f"Total gene–peak associations: {len(matches_df)}")
print(f"Unique overlapping ATAC peaks: {unique_peaks}")
print(f"Output saved to: {output_file}")
print("-" * 50)
