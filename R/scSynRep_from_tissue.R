#' @title Complete Pipeline for Synthetic Replicate Generation from Tissue
#'
#' @description A wrapper function that runs the complete scSynRep pipeline:
#'   preprocessing, bulk distribution fitting, scDesign3 modeling, and synthetic
#'   data generation.
#'
#' @param tissue_name Character; name of the tissue (matches GTEx tissue file name).
#' @param scRNA_matrix Numeric matrix of scRNA-seq counts (genes x cells) with gene names as rownames.
#' @param gtex_data_dir Character; directory containing GTEx tissue RDS files.
#' @param save_dir Character; directory to save scDesign3 model outputs.
#' @param replicate_dir Character; directory to save synthetic replicates.
#' @param Cell_label Optional DataFrame with cell annotations.
#' @param number.pc Integer; number of principal components. Default: 10.
#' @param number.gene Integer; number of highly variable genes. Default: 1500.
#' @param n_cores_bulk Integer; cores for bulk fitting. Default: 20.
#' @param use.option Integer; scDesign3 modeling option. Default: 2.
#' @param celltype_col Character; cell type column name. Default: "cell_type".
#' @param mu_formula Character; mean formula. Default: "cell_type".
#' @param n_cores_marginal Integer; cores for marginal fitting. Default: 10.
#' @param family_copula Character; copula family. Default: "nb".
#' @param parallel_para Character; parallelization method. Default: "pbmcmapply".
#' @param use.cor Integer; correlation method for bulk sampling. Default: 3.
#' @param min.eig Numeric; minimum eigenvalue. Default: 1.
#' @param number.replicate Integer; number of synthetic replicates. Default: 100.
#' @param match.option Integer; matching method. Default: 2.
#' @param scaling_factor Numeric; scaling factor. Default: 1.
#'
#' @return Invisibly returns a list with all intermediate outputs.
#'
#' @export
#' @import SummarizedExperiment
#'
#' @examples
#' \dontrun{
#' result <- scSynRep_from_tissue(
#'   tissue_name = "Liver",
#'   scRNA_matrix = my_sc_matrix,
#'   gtex_data_dir = "path/to/gtex/",
#'   save_dir = "output/model/",
#'   replicate_dir = "output/replicates/",
#'   number.replicate = 50
#' )
#' }
scSynRep_from_tissue <- function(
    tissue_name,
    scRNA_matrix,
    gtex_data_dir      = "/home/chengfeng/scRobust/data/GTEx/tissue_data/",
    save_dir,
    replicate_dir,
    Cell_label         = NULL,

    # scSynRep_prep parameters
    number.pc          = 10,
    number.gene        = 1500,

    # fit_bulk parameters
    n_cores_bulk       = 20,

    # scDesign3_fit parameters
    use.option         = 2,
    celltype_col       = "cell_type",
    mu_formula         = "cell_type",
    n_cores_marginal   = 10,
    family_copula      = "nb",
    parallel_para      = "pbmcmapply",

    # scSynRep_gen_bulk parameters
    use.cor            = 3,
    min.eig            = 1,
    number.replicate   = 100,

    # scSynRep_gen_sc parameters
    match.option       = 2,
    scaling_factor     = 1
) {

  # 0. Validate inputs
  if (is.null(tissue_name) || !is.character(tissue_name) || nchar(tissue_name) == 0) {
    stop("tissue_name must be a non-empty character string")
  }

  if (is.null(scRNA_matrix)) {
    stop("scRNA_matrix must be provided")
  }

  if (is.null(rownames(scRNA_matrix))) {
    stop("scRNA_matrix must have gene names as rownames")
  }

  # 1. Construct file path and load bulk RNA-seq data
  message("===============================================================")
  message(sprintf("Loading bulk RNA-seq data for tissue: %s", tissue_name))
  message("===============================================================")

  if (!endsWith(gtex_data_dir, "/")) {
    gtex_data_dir <- paste0(gtex_data_dir, "/")
  }

  bulk_file_path <- paste0(gtex_data_dir, tissue_name, ".RDS")

  if (!file.exists(bulk_file_path)) {
    bulk_file_path_lower <- paste0(gtex_data_dir, tolower(tissue_name), ".RDS")
    if (file.exists(bulk_file_path_lower)) {
      bulk_file_path <- bulk_file_path_lower
    } else {
      available_files <- list.files(gtex_data_dir, pattern = "\\.RDS$", ignore.case = TRUE)
      available_tissues <- gsub("\\.RDS$", "", available_files, ignore.case = TRUE)
      stop(sprintf(
        "Bulk RNA-seq file not found: %s\nAvailable tissues: %s",
        bulk_file_path,
        paste(available_tissues, collapse = ", ")
      ))
    }
  }

  bulkRNA <- readRDS(bulk_file_path)
  bulkRNA_matrix <- SummarizedExperiment::assay(bulkRNA)
  rownames(bulkRNA_matrix) <- SummarizedExperiment::rowData(bulkRNA)$Description

  message(sprintf("Loaded bulk RNA-seq matrix: %d genes x %d samples",
                  nrow(bulkRNA_matrix), ncol(bulkRNA_matrix)))

  # 2. Run scSynRep_prep
  message("\n===============================================================")
  message("Step 1: Preprocessing bulk and single-cell data")
  message("===============================================================")

  ouput1 <- scSynRep_prep(
    bulkRNA_matrix = bulkRNA_matrix,
    scRNA_matrix   = scRNA_matrix,
    number.pc      = number.pc,
    number.gene    = number.gene
  )

  message(sprintf("Selected %d highly variable genes", length(ouput1$hvg)))
  message(sprintf("Filtered bulk matrix: %d genes x %d samples",
                  nrow(ouput1$bulk), ncol(ouput1$bulk)))
  message(sprintf("Filtered sc matrix: %d genes x %d cells",
                  nrow(ouput1$sc), ncol(ouput1$sc)))

  # 3. Fit bulk distribution
  message("\n===============================================================")
  message("Step 2: Fitting bulk RNA-seq distribution")
  message("===============================================================")

  ouput2 <- fit_bulk(
    bulkRNA_matrix = ouput1$bulk,
    n_cores        = n_cores_bulk
  )

  message("Bulk distribution fitting completed")

  # 4. Fit scDesign3
  message("\n===============================================================")
  message("Step 3: Fitting scDesign3 model for single-cell data")
  message("===============================================================")

  ouput3 <- scDesign3_fit(
    scRNA_matrix     = ouput1$sc,
    top_pcs          = ouput1$pca,
    Cell_label       = Cell_label,
    save_dir         = save_dir,
    use.option       = use.option,
    celltype_col     = celltype_col,
    mu_formula       = mu_formula,
    n_cores_marginal = n_cores_marginal,
    family_copula    = family_copula,
    parallel_para    = parallel_para
  )

  message("scDesign3 fitting completed")

  # 5. Generate synthetic bulk samples
  message("\n===============================================================")
  message("Step 4: Generating synthetic bulk samples")
  message("===============================================================")

  ouput4 <- scSynRep_gen_bulk(
    bulkRNA_matrix   = ouput1$bulk,
    mu               = ouput2$mu,
    cov              = ouput2$cov,
    optimal_c        = ouput2$optimal_c,
    scRNA_matrix     = ouput1$sc,
    use.cor          = use.cor,
    number.replicate = number.replicate,
    min.eig          = min.eig
  )

  message(sprintf("Generated %d synthetic bulk replicates", number.replicate))

  # 6. Generate synthetic single-cell data
  message("\n===============================================================")
  message("Step 5: Generating synthetic single-cell replicates")
  message("===============================================================")

  ouput5 <- scSynRep_gen_sc(
    bulkRNA_matrix = ouput1$bulk,
    bulk_synth     = ouput4$sampling_bulk,
    mu             = ouput2$mu,
    d              = ouput4$d,
    optimal_c      = ouput2$optimal_c,
    scRNA_matrix   = ouput1$sc,
    scRNA_list     = ouput3,
    match.option   = match.option,
    scaling_factor = scaling_factor,
    save.dir       = replicate_dir
  )

  message(sprintf("Synthetic single-cell replicates saved to: %s", replicate_dir))

  # 7. Summary and return
  message("\n===============================================================")
  message("Pipeline completed successfully!")
  message("===============================================================")
  message(sprintf("Tissue: %s", tissue_name))
  message(sprintf("Number of replicates: %d", number.replicate))
  message(sprintf("Synthetic sc files: %s/replicate*.csv", replicate_dir))

  invisible(list(
    ouput1 = ouput1,
    ouput2 = ouput2,
    ouput3 = ouput3,
    ouput4 = ouput4,
    ouput5 = ouput5,
    tissue_name = tissue_name,
    bulkRNA_matrix = bulkRNA_matrix
  ))
}
