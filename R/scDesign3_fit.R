#' @title Fit scDesign3 Model for Single-Cell RNA-seq
#'
#' @description Fits the scDesign3 model to single-cell RNA-seq data, including
#'   marginal distributions and copula structure. Saves intermediate results
#'   to disk.
#'
#' @param scRNA_matrix Numeric matrix of scRNA-seq counts (genes x cells).
#' @param top_pcs Matrix of PCA embeddings (cells x PCs).
#' @param Cell_label Optional DataFrame with cell annotations.
#' @param save_dir Character; directory path to save model outputs.
#' @param use.option Integer; modeling option (1, 2, or 3). Default: 1.
#' @param assay_use Character; assay name. Default: "counts".
#' @param celltype_col Character; column name for cell types. Default: "cell_type".
#' @param pseudotime_col Character; column name for pseudotime. Default: NULL.
#' @param spatial_col Character; column names for spatial coordinates. Default: NULL.
#' @param other_covariates Character vector; additional covariates. Default: NULL.
#' @param corr_by Character; grouping variable for correlation. Default: "1".
#' @param predictor Character; predictor type. Default: "gene".
#' @param mu_formula Character; formula for mean model. Default: "cell_type".
#' @param sigma_formula Formula for dispersion model. Default: 1.
#' @param family_marginal Character; marginal distribution family. Default: "nb".
#' @param n_cores_marginal Integer; cores for marginal fitting. Default: 10.
#' @param usebam Logical; use bam for fitting. Default: FALSE.
#' @param parallel_marginal Character; parallelization method. Default: "pbmcmapply".
#' @param family_copula Character; copula family. Default: "nb".
#' @param copula Character; copula type. Default: "gaussian".
#' @param n_cores_copula Integer; cores for copula fitting. Default: 10.
#' @param parallel_copula Character; parallelization method. Default: "pbmcmapply".
#' @param n_cores_para Integer; cores for parameter extraction. Default: 10.
#' @param family_para Character; family for parameters. Default: "nb".
#' @param parallel_para Character; parallelization method. Default: "pbmcmapply".
#'
#' @return A list containing:
#'   \item{scRNA_sce_pc}{SingleCellExperiment object}
#'   \item{scRNA_data_pc}{Constructed data object}
#'   \item{scRNA_marginal_pc}{Fitted marginal distributions}
#'   \item{scRNA_copula_pc}{Fitted copula model}
#'   \item{scRNA_para_pc}{Extracted parameters}
#'
#' @export
#' @import scDesign3
#' @import SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#'
#' @examples
#' \dontrun{
#' sc_model <- scDesign3_fit(sc_matrix, pca_embeddings, save_dir = "output/")
#' }
scDesign3_fit <- function(
    scRNA_matrix,
    top_pcs,
    Cell_label           = NULL,
    save_dir,
    use.option           = 1,

    ## construct_data args
    assay_use            = "counts",
    celltype_col         = "cell_type",
    pseudotime_col       = NULL,
    spatial_col          = NULL,
    other_covariates     = NULL,
    corr_by              = "1",

    ## fit_marginal args
    predictor            = "gene",
    mu_formula           = "cell_type",
    sigma_formula        = 1,
    family_marginal      = "nb",
    n_cores_marginal     = 10,
    usebam               = FALSE,
    parallel_marginal    = "pbmcmapply",

    ## fit_copula args
    family_copula        = "nb",
    copula               = "gaussian",
    n_cores_copula       = 10,
    parallel_copula      = "pbmcmapply",

    ## extract_para args
    n_cores_para         = 10,
    family_para          = "nb",
    parallel_para        = "pbmcmapply"
) {
  # 1. ensure save_dir exists
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  # 2. Build SCE based on use.option
  if (use.option == 1) {
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays  = list(counts = scRNA_matrix),
      colData = Cell_label
    )
    for (i in seq_len(ncol(top_pcs))) {
      sce[[paste0("pc", i)]] <- top_pcs[, i]
    }
  }

  if (use.option == 2) {
    Cell_label <- S4Vectors::DataFrame(
      cell_type = rep(1, ncol(scRNA_matrix)),
      row.names = colnames(scRNA_matrix)
    )

    other_covariates <- paste0("pc", seq_len(ncol(top_pcs)))
    mu_formula <- paste(other_covariates, collapse = " + ")

    sce <- SingleCellExperiment::SingleCellExperiment(
      assays  = list(counts = scRNA_matrix),
      colData = Cell_label
    )
    for (i in seq_len(ncol(top_pcs))) {
      sce[[paste0("pc", i)]] <- top_pcs[, i]
    }
  }

  if (use.option == 3) {
    Cell_label <- S4Vectors::DataFrame(
      cell_type = rep(1, ncol(scRNA_matrix)),
      row.names = colnames(scRNA_matrix)
    )

    other_covariates <- paste0("pc", seq_len(ncol(top_pcs)))
    mu_formula <- paste(other_covariates, collapse = " + ")
  }

  # 4. construct data
  data_pc <- scDesign3::construct_data(
    sce               = sce,
    assay_use         = assay_use,
    celltype          = celltype_col,
    pseudotime        = pseudotime_col,
    spatial           = spatial_col,
    other_covariates  = other_covariates,
    corr_by           = corr_by
  )

  # 5. fit marginal
  message("Fitting marginal...")
  marginal_pc <- scDesign3::fit_marginal(
    data            = data_pc,
    predictor       = predictor,
    mu_formula      = mu_formula,
    sigma_formula   = sigma_formula,
    family_use      = family_marginal,
    n_cores         = n_cores_marginal,
    usebam          = usebam,
    parallelization = parallel_marginal
  )

  # 6. fit copula
  message("Fitting copula...")
  copula_pc <- scDesign3::fit_copula(
    sce             = sce,
    assay_use       = assay_use,
    marginal_list   = marginal_pc,
    family_use      = family_copula,
    copula          = copula,
    n_cores         = n_cores_copula,
    input_data      = data_pc$dat,
    parallelization = parallel_copula
  )

  # 7. extract parameters
  message("Extracting parameters...")
  para_pc <- scDesign3::extract_para(
    sce             = sce,
    marginal_list   = marginal_pc,
    n_cores         = n_cores_para,
    family_use      = family_para,
    new_covariate   = data_pc$newCovariate,
    data            = data_pc$dat,
    parallelization = parallel_para
  )

  # 8. save objects
  saveRDS(sce,         file = file.path(save_dir, "scRNA_sce_pc.rds"))
  saveRDS(data_pc,     file = file.path(save_dir, "scRNA_data_pc.rds"))
  saveRDS(marginal_pc, file = file.path(save_dir, "scRNA_marginal_pc.rds"))
  saveRDS(copula_pc,   file = file.path(save_dir, "scRNA_copula_pc.rds"))
  saveRDS(para_pc,     file = file.path(save_dir, "scRNA_para_pc.rds"))

  # 9. return list
  invisible(list(
    scRNA_sce_pc      = sce,
    scRNA_data_pc     = data_pc,
    scRNA_marginal_pc = marginal_pc,
    scRNA_copula_pc   = copula_pc,
    scRNA_para_pc     = para_pc
  ))
}
