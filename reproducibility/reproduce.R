args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(file_arg)) stop("Run this script with Rscript.")
setwd(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
workers_arg <- grep("^--workers=", args, value = TRUE)
workers <- if (length(workers_arg)) as.integer(sub("^--workers=", "", workers_arg[1])) else 1L
if (!is.finite(workers) || workers < 1L) stop("Workers must be a positive integer.")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run_r <- function(path, options = character(), directory = ".") {
  old <- getwd()
  on.exit(setwd(old))
  setwd(directory)
  status <- system2(rscript, c(shQuote(path), options))
  if (status != 0L) stop("R script failed: ", path)
}
if (!"--reuse-results" %in% args) {
  run_r("simulation/run_all_formal_chapter3.R", paste0("--workers=", workers))
  run_r("chapter4_clinical_trial_illustration/scripts/run_all_methods.R",
        c("--n-rep=10000", paste0("--workers=", workers)))
}
run_r("simulation/make_formal_chapter3_outputs.R")
run_r("simulation/prepare_plot_data.R")
run_r("chapter4_clinical_trial_illustration/scripts/rebuild_tables.R")
run_r("chapter4_clinical_trial_illustration/scripts/prepare_outputs.R")
for (section in c("simulation", "clinical_trial")) {
  source_dir <- if (section == "simulation") "simulation/formal_output/tables" else "chapter4_clinical_trial_illustration/tables"
  destination <- file.path("output", "tables", section)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(source_dir, pattern = "\\.csv$", full.names = TRUE)
  stopifnot(all(file.copy(files, destination, overwrite = TRUE)))
}
message("Reproduction complete.")
