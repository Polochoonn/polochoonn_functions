#' probeset_enrichment_analysis
#'
#' Run pathway enrichment analysis (GSEA) on a ranked DMP/probe table.
#'
#' @param DMPtable Dataframe with CpG probes (rows) and statistics (columns).
#' @param rankmetric Column to rank by ("logFC", "P.Value", "adj.P.Val"). Default "logFC".
#' @param pvalcutoffparam Nominal p-value threshold. Default 1 (no filter).
#' @param pvaltype Name of the p-value column to use. Default "P.Value".
#' @param existing_probe_ref_file Path to reference table (e.g. .rds mapping probes to pathways).
#' @param probe_type_selection_list Named list for optional probe subsetting.
#' @param genesetparam Not used, kept for compatibility.
#' @param speciesparam Not used, kept for compatibility.
#' @param seedparam Random seed. Default 12345.
#' @param customprobeset Not used, kept for compatibility.
#' @param num_cores Number of cores for parallel. Default 1.
#' @return Data.frame of GSEA results.
#' @import clusterProfiler dplyr tidyr parallel pryr
probeset_enrichment_analysis <- function(
    DMPtable, rankmetric = "logFC", pvalcutoffparam = 1, pvaltype = "P.Value",
    existing_probe_ref_file, probe_type_selection_list = NULL,
    genesetparam = "dummy", speciesparam = "Homo sapiens",
    seedparam = 12345, customprobeset = NULL, num_cores = 1
) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("clusterProfiler required")
  if (!requireNamespace("pryr", quietly = TRUE)) stop("pryr required")
  if (!requireNamespace("parallel", quietly = TRUE)) stop("parallel required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")

  message("\n[Memory] Starting analysis: ", format(pryr::mem_used(), big.mark = ","))
  
  if (!rankmetric %in% c("logFC", "P.Value", "adj.P.Val")) {
    stop("Invalid rankmetric. Choose from 'logFC', 'P.Value', or 'adj.P.Val'.")
  }

  message("Creating ranked list...")
  if (rankmetric == "logFC") {
    gseaintable <- DMPtable[DMPtable[, pvaltype] < pvalcutoffparam, "logFC", drop = FALSE]
    gsearanklist <- sort(setNames(gseaintable[, 1], rownames(gseaintable)), decreasing = TRUE)
  } else {
    gseaintable <- DMPtable[DMPtable[, pvaltype] < pvalcutoffparam, c("logFC", rankmetric), drop = FALSE]
    gseaintable[, "pstat"] <- -log10(gseaintable[, rankmetric]) * ifelse(gseaintable[, "logFC"] > 0, 1, -1)
    gsearanklist <- sort(setNames(gseaintable[, "pstat"], rownames(gseaintable)), decreasing = TRUE)
  }
  rm(gseaintable)
  gc(full = TRUE)

  message("Loading probe reference...")
  m_t2g <- as.data.frame(readRDS(existing_probe_ref_file))
  message("Reference loaded: ", nrow(m_t2g), " entries")

  if (!is.null(probe_type_selection_list)) {
    message("Filtering probes...")
    probe_check_table <- do.call(cbind, lapply(seq_along(probe_type_selection_list), function(i) {
      DMPtable[, names(probe_type_selection_list)[[i]]] %in% probe_type_selection_list[[i]]
    }))
    filtered_probes <- rownames(probe_check_table[rowSums(probe_check_table) == ncol(probe_check_table),])
    gsearanklist <- sort(gsearanklist[filtered_probes], decreasing = TRUE)
    rm(probe_check_table)
  }

  total_pathways <- length(unique(m_t2g$gs_name))
  message("Total pathways: ", total_pathways)
  
  chunk_size <- if(total_pathways > 500) {
    max(10, floor(total_pathways/50))
  } else {
    40 # default chunk size if not set
  }
  
  pathway_chunks <- split(unique(m_t2g[, 1]), 
                      ceiling(seq_along(unique(m_t2g[, 1])) / chunk_size))
  message("Using ", length(pathway_chunks), " chunks of ~", chunk_size, " pathways each")

  safe_gsea <- function(geneList, TERM2GENE, seed) {
    on.exit({
      gc(full = TRUE)
      message("[Chunk ", seed - seedparam, "] Memory after: ", format(pryr::mem_used(), big.mark = ","))
    })
    tmp <- tryCatch({
      clusterProfiler::GSEA(geneList = geneList,
           TERM2GENE = TERM2GENE,
           by = "fgsea",
           minGSSize = 100,
           maxGSSize = 100000,
           nPermSimple = 10000,  
           seed = seed,
           pvalueCutoff = 1.1,
           verbose = FALSE)
    }, error = function(e) {
      message("Error in chunk: ", e$message)
      return(NULL)
    })
    return(tmp)
  }

  message("[Memory] Pre-parallel: ", format(pryr::mem_used(), big.mark = ","))
  gsea_out_list <- parallel::mclapply(seq_along(pathway_chunks), function(i) {
    current_seed <- seedparam + i
    set.seed(current_seed)
    chunk_pathways <- pathway_chunks[[i]]
    message("[Chunk ", i, "/", length(pathway_chunks), "] Processing ", 
           length(chunk_pathways), " pathways (Seed: ", current_seed, ")")
    m_t2g_select <- m_t2g[m_t2g$gs_name %in% chunk_pathways, ]
    result <- safe_gsea(gsearanklist, m_t2g_select, current_seed)
    if(is.null(result)) {
      message("[Chunk ", i, "] Retrying with smaller chunks")
      sub_chunks <- split(chunk_pathways, ceiling(seq_along(chunk_pathways)/5))
      results <- list()
      for(subchunk in sub_chunks) {
        m_t2g_sub <- m_t2g[m_t2g$gs_name %in% subchunk, ]
        results <- c(results, list(safe_gsea(gsearanklist, m_t2g_sub, current_seed)))
      }
      return(results)
    }
    return(result)
  }, mc.cores = num_cores)

  final_results <- unlist(gsea_out_list, recursive = FALSE)
  valid_results <- final_results[!sapply(final_results, is.null)]

  message("[Memory] Post-analysis: ", format(pryr::mem_used(), big.mark = ","))

  if(length(valid_results) > 0) {
    gsea_out_df <- do.call(rbind, lapply(valid_results, as.data.frame))
    message("Analysis complete. Found ", nrow(gsea_out_df), " significant terms")
    return(gsea_out_df)
  } else {
    warning("No valid results obtained")
    return(data.frame())
  }
}
