







source("simulation/formal_chapter3_engine.R")

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(prefix, default = NULL) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", prefix, "="), "", hit[1])
}

scenario <- option_value("--scenario", "positive")
pilot <- "--pilot" %in% args
n_rep <- as.integer(option_value("--n-rep", if (pilot) "10" else FORMAL_N_REP))
chunk_size <- as.integer(option_value(
  "--chunk-size", if (pilot) as.character(n_rep) else FORMAL_CHUNK_SIZE
))
force <- "--force" %in% args
workers <- as.integer(option_value("--workers", "1"))

if (!pilot && n_rep != FORMAL_N_REP) {
  stop("Formal runs must use exactly ", FORMAL_N_REP,
       " repetitions. Use --pilot for a development run.")
}

output_root <- if (pilot) "simulation/pilot_output" else FORMAL_RESULTS_DIR
save_trace <- if (pilot) TRUE else NULL

fit <- run_formal_scenario(
  scenario = scenario,
  n_rep = n_rep,
  chunk_size = chunk_size,
  output_root = output_root,
  force = force,
  save_trace = save_trace,
  workers = workers
)

print(fit$summary, row.names = FALSE)
