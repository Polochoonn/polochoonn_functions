run_limma_dmp <- function(beta, pheno, covars, group_names, comp_name, out_dir, covariates) {
  sample_mask <- !is.na(pheno) & rowSums(is.na(covars[, covariates, drop=FALSE])) == 0
  n_samples <- sum(sample_mask)

  if (n_samples < 2 || length(unique(pheno[sample_mask])) < 2) {
    message(sprintf("[%s] Not enough samples for comparison. Skipping.", comp_name))
    return(data.frame())
  }

  pheno_used <- factor(pheno[sample_mask], levels = group_names)
  beta_used <- beta[, sample_mask, drop=FALSE]
  covars_used <- covars[sample_mask, covariates, drop=FALSE]
  design <- model.matrix(~ pheno_used + ., data = covars_used)
  message(sprintf("[%s] Running: %s vs %s | Samples: %d | Probes: %d", 
    comp_name, group_names[1], group_names[2], n_samples, nrow(beta)))

  fit <- limma::lmFit(beta_used, design)
  fit <- limma::eBayes(fit)
  res <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  res$deltaBeta <- rowMeans(beta_used[, pheno_used == group_names[1], drop = FALSE]) - 
                   rowMeans(beta_used[, pheno_used == group_names[2], drop = FALSE])
  message(sprintf("[%s] Direction: deltaBeta = mean(%s) - mean(%s)", 
  comp_name, group_names[1], group_names[2]))
  
  res$feature <- rownames(beta)[1:nrow(res)]
  res <- res %>%
    dplyr::rename(P.Value = "P.Value", adj.P.Val = "adj.P.Val") %>%
    dplyr::select(feature, logFC, AveExpr, t, P.Value, adj.P.Val, B, deltaBeta)
  rownames(res) <- rownames(beta)
  sig_count <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
  message(sprintf("[%s] Significant CpGs (adj.P.Val < 0.05): %d", comp_name, sig_count))

  out_subdir <- file.path(out_dir, comp_name)
  dir.create(out_subdir, showWarnings = FALSE, recursive = TRUE)
  write.csv(res, file = file.path(out_subdir, paste0(comp_name, "_DMP_results.csv")), row.names = TRUE)
  
  # Return only results (plotting can be handled in main script if desired)
  return(res)
}
