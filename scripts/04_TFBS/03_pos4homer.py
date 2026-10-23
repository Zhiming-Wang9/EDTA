
#!/usr/bin/env python3
"""
Step 03: Prepare target and background regions for HOMER.

Inputs:
    results/04_TFBS/{compound}_posi_CREs.tsv
    data/ATAC_GRCz11_peaks_434221.bed

Outputs:
    results/04_TFBS/{compound}_posi.homer.pos
    results/04_TFBS/ATAC_background.homer.pos

Reference genome: GRCz11
"""

import argparse
from pathlib import Path

import pandas as pd


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
    description="Prepare HOMER target and background regions."
)

parser.add_argument(
    "--compound",
    required=True,
    choices=["BET", "FLU"],
    help="Compound name: BET or FLU"
)

args = parser.parse_args()
compound = args.compound

input_file = OUTPUT_DIR / f"{compound}_posi_CREs.tsv"

target_file = OUTPUT_DIR / f"{compound}_posi.homer.pos"

background_file = OUTPUT_DIR / "ATAC_background.homer.pos"


# --------------------------------------------------
# 3. Check input files
# --------------------------------------------------

if not input_file.is_file():
    raise FileNotFoundError(
        f"CRE association table not found: {input_file}"
    )

if not ATAC_FILE.is_file():
    raise FileNotFoundError(
        f"ATAC-seq BED file not found: {ATAC_FILE}"
    )

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# 4. Prepare target ATAC-seq regions
# --------------------------------------------------

cre = pd.read_csv(input_file, sep="\t")

required_columns = [
    "chrom", "start_bed", "end_bed"
]

missing_columns = [
    col for col in required_columns
    if col not in cre.columns
]

if missing_columns:
    raise ValueError(
        f"Missing columns in CRE table: {missing_columns}"
    )

target = cre[required_columns].copy()

target["chrom"] = target["chrom"].astype(str).apply(
    lambda x: x if x.startswith("chr") else "chr" + x
)

target["start_bed"] = pd.to_numeric(
    target["start_bed"], errors="raise"
).astype(int)

target["end_bed"] = pd.to_numeric(
    target["end_bed"], errors="raise"
).astype(int)

# Retain each genomic region only once.
target = target.drop_duplicates(
    subset=["chrom", "start_bed", "end_bed"]
).reset_index(drop=True)

if target.empty:
    raise ValueError("No target ATAC-seq peaks found.")

if (
    (target["start_bed"] < 0) |
    (target["end_bed"] <= target["start_bed"])
).any():
    raise ValueError("Invalid target peak coordinates.")

target_pos = pd.DataFrame({
    "chr": target["chrom"],
    "start": target["start_bed"],
    "end": target["end_bed"],
    "name": [
        f"{compound}_peak_{i + 1}"
        for i in range(len(target))
    ],
    "score": 0,
    "strand": "+"
})

target_pos.to_csv(
    target_file,
    sep="\t",
    header=False,
    index=False
)


# --------------------------------------------------
# 5. Prepare background ATAC-seq regions
# --------------------------------------------------

# Read the first three columns of the original BED file.
# These represent chromosome, start, and end coordinates.

background = pd.read_csv(
    ATAC_FILE,
    sep="\t",
    comment="#",
    header=None,
    usecols=[0, 1, 2],
    names=["chr", "start", "end"],
    dtype={"chr": str}
)

background["chr"] = background["chr"].apply(
    lambda x: x if x.startswith("chr") else "chr" + x
)

background["start"] = pd.to_numeric(
    background["start"], errors="raise"
).astype(int)

background["end"] = pd.to_numeric(
    background["end"], errors="raise"
).astype(int)

# Remove duplicated genomic regions.
background = background.drop_duplicates(
    subset=["chr", "start", "end"]
).reset_index(drop=True)

if (
    (background["start"] < 0) |
    (background["end"] <= background["start"])
).any():
    raise ValueError("Invalid background peak coordinates.")

background["name"] = [
    f"background_peak_{i + 1}"
    for i in range(len(background))
]

background["score"] = 0
background["strand"] = "+"

background = background[
    ["chr", "start", "end", "name", "score", "strand"]
]

background.to_csv(
    background_file,
    sep="\t",
    header=False,
    index=False
)


# --------------------------------------------------
# 6. Report summary statistics
# --------------------------------------------------

print("-" * 50)
print(f"Compound: {compound}")
print(f"Input gene–peak associations: {len(cre)}")
print(f"Unique target ATAC peaks: {len(target_pos)}")
print(f"Unique background ATAC peaks: {len(background)}")
print(f"Target file: {target_file}")
print(f"Background file: {background_file}")
print("-" * 50)
