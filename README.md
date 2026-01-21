# scSynRep

[![R-CMD-check](https://github.com/JSB-UCLA/scSynRep/workflows/R-CMD-check/badge.svg)](https://github.com/JSB-UCLA/scSynRep/actions)

Generate synthetic single-cell RNA-seq samples by mapping bulk RNA-seq distribution fitting with scDesign3-based single-cell modeling.

## Overview

**scSynRep** provides a pipeline for generating realistic synthetic single-cell RNA-seq samples that capture biological variation observed in multi-sample bulk RNA-seq experiments. This is useful for:

- Prioritizing stable discoveries
- Selecting data-specific stable method

## Installation

### From GitHub

```r
# Install devtools if needed
install.packages("devtools")

# Install scSynRep
devtools::install_github("JSB-UCLA/scSynRep")
```

## Quick Start

### Basic Usage

```r
library(scSynRep)

# Load your data
bulk_matrix <- readRDS("path/to/bulk_matrix.rds")
sc_matrix <- readRDS("path/to/sc_matrix.rds")

# Step 1: Preprocess data
prep_result <- scSynRep_prep(
    bulkRNA_matrix = bulk_matrix,
    scRNA_matrix = sc_matrix,
    number.gene = 100,
    number.pc = 10
)

# Step 2: Fit bulk distribution
bulk_fit <- fit_bulk(
    bulkRNA_matrix = prep_result$bulk,
    n_cores = 10
)

# Step 3: Fit scDesign3 model
sc_model <- scDesign3_fit(
    scRNA_matrix = prep_result$sc,
    top_pcs = prep_result$pca,
    save_dir = "output/model/",
    use.option = 2,
    n_cores_marginal = 10
)

# Step 4: Generate synthetic bulk samples
synth_bulk <- scSynRep_gen_bulk(
    bulkRNA_matrix = prep_result$bulk,
    mu = bulk_fit$mu,
    cov = bulk_fit$cov,
    optimal_c = bulk_fit$optimal_c,
    scRNA_matrix = prep_result$sc,
    number.replicate = 50
)

# Step 5: Generate synthetic single-cell replicates
scSynRep_gen_sc(
    bulkRNA_matrix = prep_result$bulk,
    bulk_synth = synth_bulk$sampling_bulk,
    mu = bulk_fit$mu,
    d = synth_bulk$d,
    optimal_c = bulk_fit$optimal_c,
    scRNA_matrix = prep_result$sc,
    scRNA_list = sc_model,
    save.dir = "output/replicates/"
)
```

### Using GTEx Bulk Data (Complete Pipeline)

```r
library(scSynRep)

# Load your single-cell data
sc_matrix <- readRDS("path/to/sc_matrix.rds")

# Run complete pipeline with GTEx tissue data
result <- scSynRep_from_tissue(
    tissue_name = "Liver",
    scRNA_matrix = sc_matrix,
    gtex_data_dir = "path/to/gtex/tissue_data/",
    save_dir = "output/model/",
    replicate_dir = "output/replicates/",
    number.gene = 1500,
    number.pc = 10,
    number.replicate = 100
)
```

## Data Format Requirements

### Bulk RNA-seq Matrix
- Genes as rows, samples as columns
- Must have gene names as rownames
- Normalized expression values (e.g., TPM, FPKM)

### Single-Cell RNA-seq Matrix
- Genes as rows, cells as columns
- Must have gene names as rownames
- Raw counts (not normalized)

## Output

Synthetic replicates are saved as tab-separated CSV files:
- `replicate1.csv`, `replicate2.csv`, ..., `replicateN.csv`

Each file contains a count matrix with genes as rows and cells as columns.

## Folder Structure for Example Data

```
your_project/
├── data/
│   ├── bulk_matrix.rds        # Your bulk RNA-seq matrix
│   ├── sc_matrix.rds          # Your scRNA-seq count matrix
│   └── gtex/                   # Optional: GTEx tissue data
│       ├── Liver.RDS
│       ├── Brain.RDS
│       └── ...
├── output/
│   ├── model/                  # scDesign3 model files
│   └── replicates/             # Generated synthetic replicates
└── scripts/
    └── run_scSynRep.R
```

## License

GPL-3.0 License

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## Contact

- Chengfeng Jiang
- Project Link: https://github.com/JSB-UCLA/scSynRep
