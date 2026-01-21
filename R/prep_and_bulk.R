#' @title Filter and Preprocess Bulk and Single-Cell RNA-seq Data
#'
#' @description Filters bulk RNA-seq and scRNA-seq matrices, identifies common
#'   genes, selects highly variable genes (HVGs), and optionally computes PCA
#'   embeddings.
#'
#' @param bulkRNA_matrix Numeric matrix of bulk RNA-seq expression (genes x samples).
#'   Must have gene names as rownames.
#' @param scRNA_matrix Numeric matrix of scRNA-seq counts (genes x cells).
#'   Must have gene names as rownames.
#' @param use.pc Logical; whether to compute PCA embeddings. Default: TRUE.
#' @param number.pc Integer; number of principal components to compute. Default: 20.
#' @param number.gene Integer; number of highly variable genes to select. Default: 2000.
#' @param bulk_freq Numeric; maximum proportion of zeros allowed for bulk genes. Default: 0.2.
#' @param sc_freq Numeric; maximum proportion of zeros allowed for sc genes. Default: 0.95.
#' @param min.cells Integer; minimum number of cells a gene must be detected in. Default: 2.
#' @param select_genes Character vector; optional pre-selected genes to use instead of HVG selection.
#'
#' @return A list containing:
#'   \item{bulk}{Filtered bulk RNA-seq matrix (HVGs x samples)}
#'   \item{sc}{Filtered scRNA-seq matrix (HVGs x cells)}
#'   \item{hvg}{Character vector of highly variable gene names}
#'   \item{pca}{Matrix of PCA embeddings (cells x PCs), or NULL if use.pc=FALSE}
#'   \item{seurat}{Seurat object used for preprocessing}
#'
#' @export
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures VariableFeatures ScaleData RunPCA Embeddings
#'
#' @examples
#' \dontrun{
#' result <- scSynRep_prep(bulk_matrix, sc_matrix, number.gene = 1500)
#' }
scSynRep_prep <- function(bulkRNA_matrix,
                          scRNA_matrix,
                          use.pc       = TRUE,
                          number.pc    = 20,
                          number.gene  = 2000,
                          bulk_freq    = 0.2,
                          sc_freq      = 0.95,
                          min.cells    = 2,
                          select_genes = NULL) {
  # 0. Check inputs: both matrices must have rownames = gene names
  if (is.null(rownames(bulkRNA_matrix)) ||
      is.null(rownames(scRNA_matrix))) {
    stop("bulkRNA_matrix or scRNA_matrix must contain gene names as rownames")
  }

  # 1. Filter bulk RNA-seq genes by zero-count frequency
  min_freq <- apply(bulkRNA_matrix, 1, function(x) sum(x == min(x)) / length(x))
  indices <- which(min_freq < bulk_freq)
  keep_genes <- rownames(bulkRNA_matrix)[indices]

  min_freq <- rowMeans(scRNA_matrix == 0, na.rm = TRUE)
  indices <- which(min_freq < sc_freq)
  keep_genes2 <- rownames(scRNA_matrix)[indices]

  common_genes <- intersect(keep_genes, keep_genes2)
  if (length(common_genes) == 0) {
    stop("No common genes remain after filtering by bulk_freq")
  }
  bulkRNA_matrix_c <- bulkRNA_matrix[common_genes, , drop = FALSE]
  scRNA_matrix_c   <- scRNA_matrix[common_genes, , drop = FALSE]

  if (!is.null(select_genes)) {
    bulkRNA_matrix_c <- bulkRNA_matrix[select_genes, , drop = FALSE]
    scRNA_matrix_c   <- scRNA_matrix[select_genes, , drop = FALSE]
  }

  # 2. Create Seurat object and normalize single-cell data
  scS <- Seurat::CreateSeuratObject(
    counts       = scRNA_matrix_c,
    project      = "scS",
    min.cells    = min.cells,
    min.features = 0
  )
  scS <- Seurat::NormalizeData(scS)

  # 3. Select highly variable genes (HVGs)
  max_genes <- nrow(scS)
  if (number.gene > max_genes) {
    warning(sprintf(
      "number.gene (%d) exceeds available genes (%d); resetting to %d.",
      number.gene, max_genes, max_genes
    ))
    number.gene <- max_genes
  }
  scS <- Seurat::FindVariableFeatures(
    object           = scS,
    selection.method = "vst",
    nfeatures        = number.gene
  )

  if (is.null(select_genes)) {
    hvg <- Seurat::VariableFeatures(scS)
  } else {
    hvg <- select_genes
  }

  # 4. (Optional) Scale HVGs and run PCA
  if (use.pc) {
    scS <- Seurat::ScaleData(scS, features = hvg)
    scS <- Seurat::RunPCA(scS, features = hvg, npcs = number.pc)
    pca_embed <- Seurat::Embeddings(scS, "pca")
    n_pc_avail <- ncol(pca_embed)
    pca_mat    <- pca_embed[, seq_len(min(number.pc, n_pc_avail)), drop = FALSE]
  } else {
    pca_mat <- NULL
  }

  # 5. Subset bulk & single-cell matrices to HVGs
  bulk_out <- bulkRNA_matrix_c[hvg, , drop = FALSE]
  sc_out   <- as.matrix(scRNA_matrix_c[hvg, , drop = FALSE])

  # 6. Return a list
  list(
    bulk   = bulk_out,
    sc     = sc_out,
    hvg    = hvg,
    pca    = pca_mat,
    seurat = scS
  )
}


#' @title Fit Bulk RNA-seq Distribution
#'
#' @description Fits a multivariate normal distribution to log-transformed bulk
#'   RNA-seq data. For each gene, finds an optimal constant c that maximizes
#'   normality of log(counts + c).
#'
#' @param bulkRNA_matrix Numeric matrix of bulk RNA-seq expression (genes x samples).
#' @param n_cores Integer; number of cores for parallel processing. Default: 20.
#' @param p_val_threshold Numeric; target p-value threshold for normality test. Default: 0.5.
#' @param step Numeric; step size for searching optimal c. Default: 1.
#' @param c_value_max Numeric; maximum value of c to search. Default: 200.
#' @param normal_test Character; normality test to use ("auto", "shapiro", or "ad"). Default: "auto".
#' @param sw_limit Integer; sample size limit for Shapiro-Wilk test. Default: 200.
#'
#' @return A list containing:
#'   \item{bulk_data_counts}{Log-transformed bulk expression matrix}
#'   \item{mu}{Named vector of gene means}
#'   \item{cov}{Covariance matrix of log-transformed expression}
#'   \item{optimal_c}{Named vector of optimal c values per gene}
#'
#' @export
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom nortest ad.test
#' @importFrom stats cov shapiro.test
#'
#' @examples
#' \dontrun{
#' bulk_fit <- fit_bulk(filtered_bulk_matrix, n_cores = 10)
#' }
fit_bulk <- function(
    bulkRNA_matrix,
    n_cores         = 20,
    p_val_threshold = 0.5,
    step            = 1,
    c_value_max     = 200,
    normal_test     = c("auto"),
    sw_limit        = 200
) {
  # Define internal helper: choose c to maximize normality p-value
  find_optimal_c <- function(
    gene_counts,
    p_val_threshold,
    step,
    c_value_max,
    normal_test,
    sw_limit
  ) {
    optimal_c       <- 1
    optimal_p_value <- 0

    for (c_val in seq(1, c_value_max, by = step)) {
      log_counts <- log(gene_counts + c_val)
      n_samples  <- length(log_counts)

      # pick normality test
      if (normal_test == "shapiro" ||
          (normal_test == "auto" && n_samples <= sw_limit)) {
        p_val <- tryCatch(stats::shapiro.test(log_counts)$p.value,
                          error = function(e) NA_real_)
      } else {
        p_val <- tryCatch(nortest::ad.test(log_counts)$p.value,
                          error = function(e) NA_real_)
      }

      # update if improved
      if (!is.na(p_val) && p_val > optimal_p_value) {
        optimal_c       <- c_val
        optimal_p_value <- p_val
        if (optimal_p_value >= p_val_threshold) break
      }
    }

    list(optimal_c       = optimal_c,
         optimal_p_value = optimal_p_value)
  }

  # 1. Spin up cluster and register for foreach
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  # 2. Parallel loop over genes
  results <- foreach::foreach(
    i         = seq_len(nrow(bulkRNA_matrix)),
    .combine  = rbind,
    .packages = "nortest"
  ) %dopar% {
    find_optimal_c(
      gene_counts     = bulkRNA_matrix[i, ],
      p_val_threshold = p_val_threshold,
      step            = step,
      c_value_max     = c_value_max,
      normal_test     = normal_test,
      sw_limit        = sw_limit
    )
  }

  # 3. Tear down cluster
  parallel::stopCluster(cl)

  # 4. Assign row & column names
  rownames(results) <- rownames(bulkRNA_matrix)
  colnames(results) <- c("optimal_c", "optimal_p_value")

  # 5. Extract c-values and log-transform counts
  optimal_c        <- stats::setNames(as.numeric(results[, "optimal_c"]),
                                      rownames(results))
  bulk_data_counts <- log(bulkRNA_matrix + optimal_c)

  # 6. Compute summary statistics
  mu    <- rowMeans(bulk_data_counts)
  cov_m <- stats::cov(t(bulk_data_counts))

  # 7. Return as list
  list(
    bulk_data_counts = bulk_data_counts,
    mu               = mu,
    cov              = cov_m,
    optimal_c        = optimal_c
  )
}
