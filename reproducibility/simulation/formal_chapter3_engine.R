

source("simulation/formal_chapter3_config.R")

align_rows <- function(rows) {
  columns <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(columns, names(row))
    for (label in missing) row[[label]] <- NA
    row[columns]
  })
  do.call(rbind, rows)
}

simulate_formal_replicate <- function(scenario, replicate,
                                      save_trace = FALSE) {
  registry <- chapter3_scenarios()
  scenario_seed <- registry$seed[registry$scenario == scenario]
  config <- chapter3_config(scenario)
  gate <- NULL
  stream_seed <- scenario_seed + 100003L * as.integer(replicate)
  stream <- generate_canonical_stream(config, stream_seed)

  summaries <- vector("list", length(SEVEN_METHODS))
  traces <- if (save_trace) vector("list", length(SEVEN_METHODS)) else NULL
  for (j in seq_along(SEVEN_METHODS)) {
    method <- SEVEN_METHODS[j]
    result <- simulate_seven_method(
      method, stream, config, gate,
      assignment_seed = stream_seed + 1009L * j,
      return_trace = save_trace
    )
    if (save_trace) {
      summaries[[j]] <- result$summary
      traces[[j]] <- result$trace
      traces[[j]]$replicate <- replicate
      traces[[j]]$scenario <- scenario
    } else {
      summaries[[j]] <- result
    }
    summaries[[j]]$replicate <- replicate
    summaries[[j]]$scenario <- scenario
  }

  list(
    summary = align_rows(summaries),
    trace = if (save_trace) align_rows(traces) else NULL
  )
}

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".partial")
  saveRDS(object, temporary, compress = "xz")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) stop("Could not finalize: ", path)
  invisible(path)
}

formal_chunk_path <- function(scenario, first, last, output_root,
                              kind = c("summary", "trace")) {
  kind <- match.arg(kind)
  file.path(
    output_root, "chunks", scenario,
    sprintf("%s_%05d_%05d.rds", kind, first, last)
  )
}

run_formal_scenario <- function(scenario, n_rep = FORMAL_N_REP,
                                chunk_size = FORMAL_CHUNK_SIZE,
                                output_root = FORMAL_RESULTS_DIR,
                                force = FALSE,
                                save_trace = NULL,
                                workers = 1L) {
  validate_formal_chapter3_config()
  registry <- chapter3_scenarios()
  if (!scenario %in% registry$scenario) stop("Unknown scenario: ", scenario)
  if (is.null(save_trace)) {
    save_trace <- registry$save_monthly_trace[registry$scenario == scenario]
  }
  workers <- max(1L, as.integer(workers))
  cluster <- NULL
  if (workers > 1L) {
    cluster <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    workspace <- normalizePath(".", winslash = "/", mustWork = TRUE)
    parallel::clusterCall(cluster, setwd, workspace)
    parallel::clusterEvalQ(cluster, {
      source("simulation/formal_chapter3_engine.R")
      NULL
    })
  }

  starts <- seq.int(1L, n_rep, by = chunk_size)
  for (first in starts) {
    last <- min(n_rep, first + chunk_size - 1L)
    summary_path <- formal_chunk_path(
      scenario, first, last, output_root, "summary"
    )
    trace_path <- formal_chunk_path(
      scenario, first, last, output_root, "trace"
    )
    complete <- file.exists(summary_path) && (!save_trace || file.exists(trace_path))
    if (complete && !force) {
      message("Skipping completed chunk ", scenario, " ", first, "-", last)
      next
    }

    message("Running ", scenario, " replicates ", first, "-", last)
    fits <- if (is.null(cluster)) {
      lapply(first:last, function(replicate) {
        simulate_formal_replicate(scenario, replicate, save_trace)
      })
    } else {
      parallel::parLapply(
        cluster, first:last,
        function(replicate, scenario, save_trace) {
          simulate_formal_replicate(scenario, replicate, save_trace)
        },
        scenario = scenario, save_trace = save_trace
      )
    }
    summary_chunk <- align_rows(lapply(fits, `[[`, "summary"))
    atomic_save_rds(summary_chunk, summary_path)
    if (save_trace) {
      trace_chunk <- align_rows(lapply(fits, `[[`, "trace"))
      atomic_save_rds(trace_chunk, trace_path)
    }
  }

  aggregate_formal_scenario(
    scenario, n_rep = n_rep, output_root = output_root,
    require_trace = save_trace
  )
}

summarise_trace <- function(trace) {
  keys <- unique(trace[c("method", "month")])
  keys <- keys[order(match(keys$method, SEVEN_METHODS), keys$month), ]
  metrics <- c(
    "randomization_probability", "target_1", "evidence_z",
    "cumulative_allocation_1", "available_z", "available_y", "enrolled"
  )
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    x <- trace[
      trace$method == keys$method[i] & trace$month == keys$month[i],
      , drop = FALSE
    ]
    answer <- keys[i, , drop = FALSE]
    for (metric in metrics) {
      value <- x[[metric]]
      value <- value[is.finite(value)]
      answer[[paste0(metric, "_mean")]] <- if (length(value)) mean(value) else NA
      answer[[paste0(metric, "_q10")]] <- if (length(value)) {
        unname(quantile(value, 0.10, names = FALSE))
      } else NA
      answer[[paste0(metric, "_q90")]] <- if (length(value)) {
        unname(quantile(value, 0.90, names = FALSE))
      } else NA
    }
    answer
  })
  do.call(rbind, rows)
}

aggregate_formal_scenario <- function(scenario, n_rep = FORMAL_N_REP,
                                      output_root = FORMAL_RESULTS_DIR,
                                      require_trace = FALSE) {
  chunk_dir <- file.path(output_root, "chunks", scenario)
  summary_files <- list.files(
    chunk_dir, pattern = "^summary_[0-9]{5}_[0-9]{5}\\.rds$",
    full.names = TRUE
  )
  if (!length(summary_files)) stop("No chunks found for ", scenario)
  results <- align_rows(lapply(summary_files, readRDS))
  results <- results[results$replicate <= n_rep, ]
  expected_rows <- n_rep * length(SEVEN_METHODS)
  if (nrow(results) != expected_rows ||
      length(unique(results$replicate)) != n_rep ||
      anyDuplicated(results[c("replicate", "method")])) {
    stop("Incomplete or duplicated formal results for ", scenario)
  }
  results <- results[order(results$replicate,
                           match(results$method, SEVEN_METHODS)), ]

  scenario_dir <- file.path(output_root, scenario)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  atomic_save_rds(results, file.path(scenario_dir, "replications.rds"))
  write.csv(results, file.path(scenario_dir, "replications.csv"), row.names = FALSE)

  config <- chapter3_config(scenario)
  summary <- summarise_seven_method_comparison(results, config)
  summary$scenario <- scenario
  write.csv(summary, file.path(scenario_dir, "summary.csv"), row.names = FALSE)

  if (require_trace) {
    trace_files <- list.files(
      chunk_dir, pattern = "^trace_[0-9]{5}_[0-9]{5}\\.rds$",
      full.names = TRUE
    )
    if (!length(trace_files)) stop("Monthly trace is missing for ", scenario)
    trace <- align_rows(lapply(trace_files, readRDS))
    trace <- trace[trace$replicate <= n_rep, ]
    expected_trace_rows <- n_rep * length(SEVEN_METHODS) * config$n_months
    if (nrow(trace) != expected_trace_rows) {
      stop("Incomplete monthly trace for ", scenario)
    }
    trace_summary <- summarise_trace(trace)
    write.csv(
      trace_summary, file.path(scenario_dir, "monthly_trace_summary.csv"),
      row.names = FALSE
    )
  }
  invisible(list(results = results, summary = summary))
}
