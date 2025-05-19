source("https://raw.githubusercontent.com/Polochoonn/polochoonn_functions/main/psea/psea_enrichment_function.R")
needed_pkgs <- c("clusterProfiler", "pryr", "parallel", "dplyr", "tidyr")
missing <- needed_pkgs[!sapply(needed_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) stop("Install required packages: ", paste(missing, collapse = ", "))

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) stop("Usage: Rscript run_psea.R <dmp_result_file> <output_folder> <ref_file>")
dmp_file <- args[1]
output_folder <- args[2]
ref_file <- args[3]
if(!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

message("Run started: ", Sys.time())
message("Input: ", dmp_file)
message("Output dir: ", output_folder)
message("Reference file: ", ref_file)

start_time <- Sys.time()
tryCatch({
  dmp_df <- read.csv(dmp_file, row.names = 1)
  res <- probeset_enrichment_analysis(
    DMPtable = dmp_df,
    existing_probe_ref_file = ref_file,
    num_cores = 4
  )
  output_file <- file.path(output_folder, paste0("PSEA_results_GO_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
  write.csv(res, output_file, row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_folder, "sessionInfo.txt"))
}, error = function(e) {
  message("Error: ", e$message)
  quit(status = 1)
})
end_time <- Sys.time()
message("Run finished: ", end_time, " (Elapsed: ", round(difftime(end_time, start_time, units = "mins"), 2), " min)")
