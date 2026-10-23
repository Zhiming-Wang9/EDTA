# EDTA: Effect-Dependent Transcriptome Analysis

## Overview

This repository organizes the analysis code used in the EDTA (Effect-Dependent Transcriptome Analysis) study to support computational reproducibility. The currently documented module identifies ATAC-seq regions overlapping selected zebrafish genes and performs transcription factor binding motif enrichment using HOMER.

**Status:** Only Module 04 (CRE association and TFBS motif enrichment) is documented here. Other modules will be added after verification.

## Repository structure

```text
EDTA/
├── data/
│   ├── BET_posi.txt
│   ├── FLU_posi.txt
│   ├── Danio_rerio.GRCz11.115.gtf
│   └── ATAC_GRCz11_peaks_434221.bed
├── scripts/
│   ├── 01_DEG/
│   ├── 02_Correlation/
│   ├── 03_KEGG/
│   ├── 04_TFBS/
│   │   ├── 01_gene_info_from_gtf_by_gene_name.py
│   │   ├── 02_gene2CREs.py
│   │   ├── 03_pos4homer.py
│   │   └── 04_run_homer.slurm
│   └── 05_Celloracle/
├── results/
├── environment/
│   ├── meme_environment.yml
│   └── meme_environment_full.yml
└── README.md
```

The `data/` directory lists the required local inputs. Availability, provenance, and redistribution permissions for reference data should be confirmed before distributing those files through GitHub.

## Module 04: CRE association and TFBS motif enrichment

### Input data

| File | Description |
| --- | --- |
| `BET_posi.txt` | Selected betamethasone (BET) gene identifiers, one per line |
| `FLU_posi.txt` | Selected fludrocortisone (FLU) gene identifiers, one per line |
| `Danio_rerio.GRCz11.115.gtf` | Zebrafish gene annotations for the GRCz11 assembly |
| `ATAC_GRCz11_peaks_434221.bed` | ATAC-seq accessible chromatin regions mapped to GRCz11 |

### Software environment

The analyses were run on a Linux high-performance computing (HPC) cluster using Conda-managed software. The `environment/` directory contains two exports of the analysis environment:

- `meme_environment.yml`: Conda environment specification exported from explicitly requested packages (`--from-history`).
- `meme_environment_full.yml`: Full export of the installed Conda environment for recording software dependencies and versions.

The analysis requires Python 3 with `pandas`, `pybedtools`, and `tqdm`; BEDTools; HOMER with the `danRer11` genome installed; and SLURM for the provided batch script.

To create an environment from the portable specification, run:

```bash
conda env create -f environment/meme_environment.yml
conda activate meme
```

The environment name and any cluster-specific paths may need adjustment. A `--from-history` export may omit packages installed with `pip` or otherwise not explicitly recorded by Conda. Installation of HOMER and its genome resources should therefore be confirmed separately. The full environment export is an archival record rather than a guaranteed portable installation file.

### Analysis workflow

**Step 01 — Gene annotation:** Match input identifiers to `gene_id` or `gene_name` in the GRCz11 GTF and export gene coordinates. Coordinates remain 1-based and inclusive in this output.

**Step 02 — Gene–ATAC associations:** Convert gene positions to 0-based, half-open BED coordinates and use BEDTools `intersect -wa -wb` to retain ATAC peaks overlapping each gene body. These are gene-body overlaps, **not** specifically TSS-centered or promoter-only peaks.

**Step 03 — HOMER input preparation:** Export each unique associated ATAC peak as a six-column position file for HOMER. The script also generates `ATAC_background.homer.pos` from all unique ATAC peaks for recordkeeping or alternative analyses; **this background file is not used in the original-parameter HOMER run described below**.

**Step 04 — Motif enrichment:** Run HOMER against the `danRer11` reference genome with `-size given -p 8 -mask`, without a user-specified `-bg` option. Thus, HOMER uses its internally generated/default background. The `-size given` option retains the supplied genomic interval lengths.

### Run commands

Run the following commands from the repository root:

```bash
cd ~/projects/EDTA

# Step 01: Gene annotations
python scripts/04_TFBS/01_gene_info_from_gtf_by_gene_name.py --compound BET
python scripts/04_TFBS/01_gene_info_from_gtf_by_gene_name.py --compound FLU

# Step 02: Gene–ATAC overlaps
python scripts/04_TFBS/02_gene2CREs.py --compound BET
python scripts/04_TFBS/02_gene2CREs.py --compound FLU

# Step 03: Target positions (and optional ATAC background file)
python scripts/04_TFBS/03_pos4homer.py --compound BET
python scripts/04_TFBS/03_pos4homer.py --compound FLU

# Step 04: HOMER analysis on SLURM
mkdir -p logs
sbatch --export=ALL,COMPOUND=BET scripts/04_TFBS/04_run_homer.slurm
sbatch --export=ALL,COMPOUND=FLU scripts/04_TFBS/04_run_homer.slurm
```

**Cluster note:** The SLURM partition and Conda activation path in `04_run_homer.slurm` are specific to the originating cluster and must be adjusted on other systems. The script must begin with `#!/bin/bash`.

### Expected output

All outputs are written under `results/04_TFBS/`.

| Output | Description |
| --- | --- |
| `BET_posit_gene_table.tsv` / `FLU_posit_gene_table.tsv` | Matched gene annotations |
| `BET_unmatched_genes.txt` / `FLU_unmatched_genes.txt` | Unmatched input identifiers |
| `BET_posi_CREs.tsv` / `FLU_posi_CREs.tsv` | Gene–ATAC peak overlap records |
| `BET_posi.homer.pos` / `FLU_posi.homer.pos` | Unique target ATAC regions for HOMER |
| `ATAC_background.homer.pos` | Optional full ATAC peak background file; not used in the original-parameter HOMER run |
| `BET_posi_motif/` / `FLU_posi_motif/` | HOMER motif enrichment results |

### Reproducibility considerations

The HOMER command used for the current reproduction follows the archived BET run parameters (`-size given -p 8 -mask`, without `-bg`). The revised preprocessing scripts nevertheless differ from the original scripts in two respects: (i) GTF gene starts are converted correctly to BED coordinates during overlap analysis; and (ii) duplicate ATAC peak intervals are removed when preparing HOMER targets, with peak-specific identifiers assigned in place of gene identifiers. These changes can affect the selected input regions and the resulting motif enrichment statistics. Accordingly, numerical equivalence to the manuscript results should be assessed against the original input files and archived output before claiming exact reproduction. Matching command-line parameters alone does not establish equivalent results.

## Citation

Citation details will be added after the manuscript metadata are finalized.

## License

A license has not yet been selected.
