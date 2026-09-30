

source("simulation/formal_chapter3_engine.R")

args <- commandArgs(trailingOnly = TRUE)
worker_arg <- grep("^--workers=", args, value = TRUE)
workers <- if (length(worker_arg)) {
  as.integer(sub("^--workers=", "", worker_arg[1]))
} else 1L

for (scenario in chapter3_scenarios()$scenario) {
  run_formal_scenario(
    scenario = scenario,
    n_rep = FORMAL_N_REP,
    chunk_size = FORMAL_CHUNK_SIZE,
    output_root = FORMAL_RESULTS_DIR,
    workers = workers
  )
}
