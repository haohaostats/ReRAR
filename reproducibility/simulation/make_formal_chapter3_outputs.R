


source("simulation/formal_chapter3_config.R")

pilot_output <- identical(Sys.getenv("CHAPTER3_PILOT_OUTPUT"), "1")
input_root <- if (pilot_output) "simulation/pilot_output" else FORMAL_RESULTS_DIR
expected_n <- if (pilot_output) {
  as.integer(Sys.getenv("CHAPTER3_PILOT_N", "2"))
} else FORMAL_N_REP
table_dir <- if (pilot_output) {
  "simulation/pilot_output/formal_output/tables"
} else FORMAL_TABLE_DIR
figure_dir <- if (pilot_output) {
  "simulation/pilot_output/formal_output/figures"
} else FORMAL_FIGURE_DIR

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_formal_results <- function(scenario) {
  path <- file.path(input_root, scenario, "replications.rds")
  if (!file.exists(path)) stop("Missing formal result: ", path)
  x <- readRDS(path)
  if (length(unique(x$replicate)) != expected_n ||
      nrow(x) != expected_n * length(SEVEN_METHODS)) {
    stop("Scenario ", scenario, " is not a complete ", expected_n,
         "-replication formal run.")
  }
  x
}

read_formal_summary <- function(scenario) {
  path <- file.path(input_root, scenario, "summary.csv")
  if (!file.exists(path)) stop("Missing formal summary: ", path)
  x <- read.csv(path, check.names = FALSE)
  if (!setequal(x$method, SEVEN_METHODS)) {
    stop("Formal summary does not contain all seven methods: ", scenario)
  }
  x$method <- factor(x$method, levels = SEVEN_METHODS)
  x[order(x$method), ]
}

all_results <- setNames(
  lapply(chapter3_scenarios()$scenario, read_formal_results),
  chapter3_scenarios()$scenario
)
all_summaries <- setNames(
  lapply(chapter3_scenarios()$scenario, read_formal_summary),
  chapter3_scenarios()$scenario
)

write_table <- function(x, stem, digits = 4L) {
  write.csv(x, file.path(table_dir, paste0(stem, ".csv")), row.names = FALSE, na = "")
}

base <- chapter3_config("positive")
table_1 <- data.frame(
  parameter = c(
    "Total sample size", "Accrual", "Control final response",
    "Experimental final response", "Investigator assessment delay",
    "Additional central-review delay", "Control sensitivity",
    "Control specificity", "Experimental sensitivity",
    "Experimental specificity", "Allocation envelope"
  ),
  value = c(
    as.character(base$n),
    paste(base$patients_per_month, "patients/month for", base$n_months, "months"),
    format(base$p_y["0"], digits = 3),
    format(base$p_y["1"], digits = 3),
    paste(base$investigator_delay_months, "months"),
    paste("median", base$central_delay_median_months,
          "months; 90th percentile", base$central_delay_q90_months, "months"),
    format(base$sensitivity["0"], digits = 4),
    format(base$specificity["0"], digits = 4),
    format(base$sensitivity["1"], digits = 4),
    format(base$specificity["1"], digits = 4),
    paste0("[", base$power_envelope_lower, ", ",
           base$power_envelope_upper, "]")
  ),
  rationale = c(
    "Clinically anchored trial size", "Fixed monthly accrual",
    "Primary data-generating mechanism", "Primary data-generating mechanism",
    "Investigator assessment schedule", "Clinically anchored review timing",
    rep("ZUMA-7 revision matrix", 4), "Pre-specified power constraint"
  ),
  stringsAsFactors = FALSE
)
write_table(table_1, "table1", digits = 4)

failure_reduction_vs_cr <- function(results) {
  cr <- results[results$method == "CR", c("replicate", "expected_failures")]
  names(cr)[2] <- "cr_expected_failures"
  rows <- lapply(SEVEN_METHODS, function(method) {
    mm <- results[results$method == method,
                  c("replicate", "expected_failures")]
    merged <- merge(mm, cr, by = "replicate")
    d <- merged$cr_expected_failures - merged$expected_failures
    data.frame(
      method = method,
      failure_reduction_vs_cr = mean(d),
      mcse = sd(d) / sqrt(length(d)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


positive_summary <- all_summaries$positive
positive_reduction <- failure_reduction_vs_cr(all_results$positive)
table_2 <- merge(
  positive_summary[c(
    "method", "allocation_1", "expected_failures",
    "contrast_stabilized_power", "contrast_stabilized_coverage"
  )],
  positive_reduction[c("method", "failure_reduction_vs_cr")],
  by = "method", sort = FALSE
)
table_2 <- table_2[match(SEVEN_METHODS, table_2$method), ]
names(table_2) <- c(
  "Method", "Allocation to experimental", "Expected failures",
  "Power", "Coverage", "Failures avoided versus CR"
)
write_table(table_2, "table2")


null_summary <- all_summaries$null
table_3 <- null_summary[c(
  "method", "allocation_1", "contrast_stabilized_power",
  "contrast_stabilized_bias", "contrast_stabilized_empirical_se",
  "contrast_stabilized_mean_se", "contrast_stabilized_coverage"
)]
names(table_3) <- c(
  "Method", "Allocation to experimental", "Type I error", "Bias",
  "Empirical SE", "Mean estimated SE", "Coverage"
)
write_table(table_3, "table3")


reversal_summary <- all_summaries$reversal
wrong_direction <- vapply(SEVEN_METHODS, function(method) {
  x <- all_results$reversal
  mean(x$allocation_1[x$method == method] > 0.5)
}, numeric(1))
table_4 <- data.frame(
  Method = SEVEN_METHODS,
  `Allocation to truly better arm` =
    1 - reversal_summary$allocation_1,
  `Probability of majority allocation to inferior arm` = wrong_direction,
  `Expected failures` = reversal_summary$expected_failures,
  Coverage = reversal_summary$contrast_stabilized_coverage,
  check.names = FALSE
)
write_table(table_4, "table4")


sensitivity_ids <- setdiff(chapter3_scenarios()$scenario,
                           c("null", "reversal"))
table_5 <- do.call(rbind, lapply(sensitivity_ids, function(scenario) {
  x <- all_summaries[[scenario]]
  data.frame(
    scenario = scenario, method = as.character(x$method),
    allocation_1 = x$allocation_1,
    expected_failures = x$expected_failures,
    power = x$contrast_stabilized_power,
    coverage = x$contrast_stabilized_coverage,
    stringsAsFactors = FALSE
  )
}))
write_table(table_5, "table5")

table_6 <- do.call(rbind, lapply(names(all_results), function(scenario) {
  x <- all_summaries[[scenario]]
  data.frame(
    scenario = scenario, method = as.character(x$method),
    expected_failures = x$expected_failures,
    realized_failures = x$realized_failures,
    stringsAsFactors = FALSE
  )
}))
write_table(table_6, "table6")

table_7 <- do.call(rbind, lapply(names(all_summaries), function(scenario) {
  x <- all_summaries[[scenario]]
  data.frame(
    scenario = scenario, method = as.character(x$method),
    bias = x$contrast_stabilized_bias,
    empirical_se = x$contrast_stabilized_empirical_se,
    mean_estimated_se = x$contrast_stabilized_mean_se,
    coverage = x$contrast_stabilized_coverage,
    rejection_probability = x$contrast_stabilized_power,
    stringsAsFactors = FALSE
  )
}))
write_table(table_7, "table7")

write.csv(chapter3_scenario_parameters(), file.path(table_dir, "scenario_registry.csv"), row.names = FALSE)
