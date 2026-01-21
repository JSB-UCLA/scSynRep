#' @title Generate Synthetic Bulk RNA-seq Samples
#'
#' @description Generates synthetic bulk RNA-seq samples by sampling from a
#'   truncated multivariate normal distribution fitted to the original bulk data.
#'
#' @param bulkRNA_matrix Original bulk RNA-seq matrix (genes x samples).
#' @param mu Named vector of gene means from fit_bulk().
#' @param cov Covariance matrix from fit_bulk().
#' @param optimal_c Named vector of optimal c values from fit_bulk().
#' @param scRNA_matrix Single-cell RNA-seq matrix (genes x cells).
#' @param use.cor Integer; correlation method (1=truncated MVN, 2=independent truncated, 3=MVN). Default: 1.
#' @param min.eig Numeric; minimum eigenvalue for positive definiteness. Default: 1e-2.
#' @param number.replicate Integer; number of synthetic samples to generate. Default: 10.
#'
#' @return A list containing:
#'   \item{d}{Shift factor between pseudo-bulk and real bulk}
#'   \item{sampling_bulk}{Matrix of synthetic bulk samples (genes x replicates)}
#'
#' @export
#' @importFrom tmvtnorm rtmvnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom corpcor make.positive.definite
#' @importFrom MASS mvrnorm
#'
#' @examples
#' \dontrun{
#' synth_bulk <- scSynRep_gen_bulk(bulk_matrix, mu, cov, optimal_c, sc_matrix)
#' }
scSynRep_gen_bulk <- function(bulkRNA_matrix, mu, cov, optimal_c, scRNA_matrix,
                              use.cor = 1, min.eig = 1e-2, number.replicate = 10) {
  pseudo_bulk <- rowSums(scRNA_matrix)
  names(pseudo_bulk) <- rownames(scRNA_matrix)
  pseudo_bulk <- pseudo_bulk * sum(bulkRNA_matrix) / ncol(bulkRNA_matrix) / sum(pseudo_bulk)
  log_pseudo_bulk <- log(pseudo_bulk + optimal_c)
  d <- mu - log_pseudo_bulk
  lower_bound <- d + log(optimal_c)
  upper_bound <- rep(Inf, length(lower_bound))

  if (use.cor == 1) {
    cov_p <- corpcor::make.positive.definite(cov, tol = min.eig)
    sampling_bulk <- t(tmvtnorm::rtmvnorm(
      n = number.replicate, mean = mu, sigma = cov_p,
      lower = lower_bound, upper = upper_bound,
      algorithm = "gibbs", start = mu
    ))
  } else if (use.cor == 2) {
    cov_p <- diag(cov)
    samples_list <- lapply(seq_along(mu), function(i) {
      truncnorm::rtruncnorm(
        n    = number.replicate,
        a    = lower_bound[i],
        b    = upper_bound[i],
        mean = mu[i],
        sd   = sqrt(cov_p[i])
      )
    })
    samples_mat <- do.call(cbind, samples_list)
    sampling_bulk <- t(samples_mat)
  } else if (use.cor == 3) {
    sampling_bulk <- replicate(
      number.replicate,
      MASS::mvrnorm(n = 1, mu = mu, Sigma = cov)
    )
  }

  rownames(sampling_bulk) <- names(mu)
  colnames(sampling_bulk) <- paste0("rep", seq_len(number.replicate))

  list(d = d, sampling_bulk = sampling_bulk)
}


#' @title Generate Synthetic Single-Cell RNA-seq Replicates
#'
#' @description Generates synthetic single-cell RNA-seq replicates by adjusting
#'   scDesign3 parameters based on synthetic bulk samples.
#'
#' @param bulkRNA_matrix Original bulk RNA-seq matrix (genes x samples).
#' @param bulk_synth Matrix of synthetic bulk samples from scSynRep_gen_bulk().
#' @param mu Named vector of gene means from fit_bulk().
#' @param d Shift factor from scSynRep_gen_bulk().
#' @param optimal_c Named vector of optimal c values from fit_bulk().
#' @param scRNA_matrix Original single-cell RNA-seq matrix (genes x cells).
#' @param scRNA_list List of scDesign3 fitted objects from scDesign3_fit().
#' @param match.option Integer; matching method (1=precise, 2=alternative). Default: 1.
#' @param scaling_factor Numeric; scaling power for per-gene mapping. Default: 1.
#' @param use.pc Integer; whether to use PC-based correction. Default: 1.
#' @param n_cores Integer; number of cores for parallel generation. Default: 1.
#' @param sc_quantile Numeric; quantile threshold for outliers. Default: 0.995.
#' @param save.dir Character; directory to save synthetic replicates.
#'
#' @return Invisibly returns NULL. Synthetic replicates are saved as CSV files.
#'
#' @export
#' @import scDesign3
#' @import BiocParallel
#' @importFrom stats IQR median quantile
#'
#' @examples
#' \dontrun{
#' scSynRep_gen_sc(bulk_matrix, synth_bulk$sampling_bulk, mu, synth_bulk$d,
#'                 optimal_c, sc_matrix, sc_model, save.dir = "replicates/")
#' }
scSynRep_gen_sc <- function(bulkRNA_matrix, bulk_synth, mu, d, optimal_c,
                            scRNA_matrix, scRNA_list, match.option = 1,
                            scaling_factor = 1, use.pc = 1, n_cores = 1,
                            sc_quantile = 0.995, save.dir) {
  if (!dir.exists(save.dir)) dir.create(save.dir, recursive = TRUE)

  pseudo_bulk <- rowSums(scRNA_matrix)
  e_v <- stats::quantile(pseudo_bulk, 0.75, na.rm = TRUE) + 1.5 * stats::IQR(pseudo_bulk, na.rm = TRUE)
  names(pseudo_bulk) <- rownames(scRNA_matrix)
  alpha <- sum(bulkRNA_matrix) / ncol(bulkRNA_matrix) / sum(pseudo_bulk)

  # Create mapping vector
  per_gene_mapping <- list()

  if (match.option == 1) {
    for (i in seq_len(ncol(bulk_synth))) {
      per_gene_mapping[[i]] <- (exp(bulk_synth[, i] - d) - optimal_c) / (exp(mu - d) - optimal_c)
      e_v_m <- stats::quantile(per_gene_mapping[[i]], 0.75, na.rm = TRUE) +
               1.5 * stats::IQR(per_gene_mapping[[i]], na.rm = TRUE)
      per_gene_mapping[[i]][per_gene_mapping[[i]] > e_v_m] <- e_v_m
      per_gene_mapping[[i]] <- per_gene_mapping[[i]] ^ scaling_factor
    }
  }

  if (match.option == 2) {
    for (i in seq_len(ncol(bulk_synth))) {
      v <- bulk_synth[, i]
      pseudo_bulk_hat <- (exp(v - d) - optimal_c) / alpha
      pseudo_bulk_hat[pseudo_bulk_hat < 0] <- 0
      pseudo_bulk_hat[pseudo_bulk_hat > e_v] <- e_v
      per_gene_mapping[[i]] <- pseudo_bulk_hat / pseudo_bulk
      e_v_m <- stats::quantile(per_gene_mapping[[i]], 0.75, na.rm = TRUE) +
               1.5 * stats::IQR(per_gene_mapping[[i]], na.rm = TRUE)
      per_gene_mapping[[i]][per_gene_mapping[[i]] > e_v_m] <- e_v_m
      per_gene_mapping[[i]] <- per_gene_mapping[[i]] ^ scaling_factor
    }
  }

  # Generate new replicates
  Generate_multiple_counts <- lapply(seq_len(ncol(bulk_synth)), function(i) {
    set.seed(i)
    mean_mat <- scRNA_list$scRNA_para_pc$mean_mat

    per_gene_mapping[[i]][per_gene_mapping[[i]] == 0] <- 1e-6
    mean_mat_weighted <- sweep(mean_mat, 2, per_gene_mapping[[i]], `*`)

    colnames(mean_mat_weighted) <- colnames(scRNA_list$scRNA_para_pc$mean_mat)

    newcount <- scDesign3::simu_new(
      sce = scRNA_list$scRNA_sce_pc,
      mean_mat = mean_mat_weighted,
      sigma_mat = scRNA_list$scRNA_para_pc$sigma_mat,
      zero_mat = scRNA_list$scRNA_para_pc$zero_mat,
      quantile_mat = NULL,
      copula_list = scRNA_list$scRNA_copula_pc$copula_list,
      n_cores = n_cores,
      family_use = "nb",
      input_data = scRNA_list$scRNA_data_pc$dat,
      new_covariate = scRNA_list$scRNA_data_pc$newCovariate,
      parallelization = "bpmcmapply",
      BPPARAM = BiocParallel::MulticoreParam(),
      important_feature = scRNA_list$scRNA_copula_pc$important_feature,
      filtered_gene = scRNA_list$scRNA_data_pc$filtered_gene
    )

    if (use.pc) {
      rs_new <- rowSums(newcount, na.rm = TRUE)
      rs_sc  <- rowSums(scRNA_matrix, na.rm = TRUE)

      fc <- rs_new / (rs_sc + 1e-6)
      fc_reverse <- rs_sc / (rs_new + 1e-6)
      fc_reverse <- fc_reverse[rs_new != 0]

      genes_fc2_high <- names(fc)[fc > 3 * max(fc_reverse) & rs_new > stats::median(rs_new)]
      newcount[genes_fc2_high, ] <- newcount[genes_fc2_high, ] * fc_reverse[genes_fc2_high]
    }

    filename <- file.path(save.dir, paste0("replicate", i, ".csv"))
    utils::write.table(newcount, filename, sep = "\t", row.names = TRUE, col.names = TRUE)

    rm(newcount)
    rm(mean_mat)
    gc()

    invisible(NULL)
  })

  invisible(NULL)
}
