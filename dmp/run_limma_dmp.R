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
  
  #Maybe in the future do that for avoid coefficient =2 
  #group_coef <- which(colnames(design) == "pheno_used" | colnames(design) == paste0("pheno_used", group_names[1]))
  #fit <- limma::lmFit(beta_used, design)
  #fit <- limma::eBayes(fit)
  #res <- limma::topTable(fit, coef = group_coef, number = Inf, sort.by = "none")
  


  
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
  

  # Plotting as before (reuses your plotting code)
  if (nrow(res) > 0) {
    volcano_plot <- EnhancedVolcano(
      res, lab = rownames(res),
      x = "deltaBeta", y = "adj.P.Val",
      title = paste0(comp_name, " Differential Methylation"),
      pCutoff = 0.1, FCcutoff = 0, xlim = c(-0.3, 0.3), ylim = c(0, 10),
      col = c("grey30", "forestgreen", "royalblue", "red2"), legendPosition = "right"
    )
    ggsave(file.path(out_subdir, paste0(comp_name, "_volcano.png")), plot = volcano_plot, width = 8, height = 6)
    hist_adj <- ggplot(res, aes(x = adj.P.Val)) +
      geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
      theme_minimal() + ggtitle(paste(comp_name, "Distribution of adj.P.Val")) + xlab("adj.P.Val") + ylab("Frequency")
    ggsave(file.path(out_subdir, paste0(comp_name, "_adj_pvalue_histogram.png")), hist_adj, width = 8, height = 6)
    hist_raw <- ggplot(res, aes(x = P.Value)) +
      geom_histogram(binwidth = 0.01, fill = "lightcoral", color = "black") +
      theme_minimal() + ggtitle(paste(comp_name, "Histogram of Raw P-values")) + xlab("P.Value") + ylab("Frequency")
    ggsave(file.path(out_subdir, paste0(comp_name, "_raw_pvalue_histogram.png")), hist_raw, width = 8, height = 6)
    density_plot <- ggplot(res, aes(x = deltaBeta)) +
      geom_density(fill = "lightblue") + theme_minimal() + ggtitle(paste(comp_name, "Density of deltaBeta")) + xlab("deltaBeta") + ylab("Density")
    ggsave(file.path(out_subdir, paste0(comp_name, "_deltaBeta_density.png")), density_plot, width = 8, height = 6)
    # QQ plot
    res$logP <- -log10(res$P.Value)
    qq_plot <- ggplot(res, aes(sample = logP)) +
      stat_qq(distribution = qunif) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      theme_minimal() + ggtitle(paste(comp_name, "QQ Plot of -log10 Raw P-values")) +
      xlab("Theoretical Quantiles (-log10)") + ylab("Sample Quantiles (-log10)")
    ggsave(file.path(out_subdir, paste0(comp_name, "_qqplot_raw_pvalues.png")), qq_plot, width = 8, height = 6)
  }
  return(res)
}
