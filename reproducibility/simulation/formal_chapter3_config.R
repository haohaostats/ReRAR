

source("simulation/seven_method_comparison.R")

FORMAL_N_REP <- 10000L
FORMAL_CHUNK_SIZE <- 100L
FORMAL_RESULTS_DIR <- "simulation/formal_results"
FORMAL_TABLE_DIR <- "simulation/formal_output/tables"
FORMAL_FIGURE_DIR <- "simulation/formal_output/figures"

chapter3_scenarios <- function() {
  data.frame(
    scenario = c(
      "positive", "null", "reversal",
      "delay_short", "delay_long",
      "transport_moderate", "transport_severe",
      "surrogate_weak", "surrogate_strong"
    ),
    family = c(
      "confirmatory", "confirmatory", "confirmatory",
      "delay", "delay", "transport", "transport",
      "surrogate", "surrogate"
    ),
    level = c(
      "anchor", "null", "reversal",
      "short", "long", "moderate", "severe", "weak", "strong"
    ),
    seed = c(
      31000001L, 32000003L, 33000007L,
      34000009L, 35000011L,
      36000013L, 37000017L,
      38000021L, 39000023L
    ),
    save_monthly_trace = c(TRUE, rep(FALSE, 8L)),
    stringsAsFactors = FALSE
  )
}

mix_arm_revision <- function(x, amount) {
  stopifnot(amount >= 0, amount <= 1, length(x) == 2L)
  answer <- (1 - amount) * x + amount * rev(x)
  names(answer) <- names(x)
  answer
}

weaken_binary_accuracy <- function(x) {
  answer <- 0.5 + 0.5 * (x - 0.5)
  names(answer) <- names(x)
  answer
}

strengthen_binary_accuracy <- function(x) {
  answer <- x + 0.5 * (1 - x)
  names(answer) <- names(x)
  answer
}

validation_counts_at_accuracy <- function(config) {
  validation_sizes <- data.frame(
    arm = 0:1,
    se_total = c(90, 150),
    sp_total = c(89, 30)
  )
  data.frame(
    arm = validation_sizes$arm,
    se_success = validation_sizes$se_total * unname(config$sensitivity),
    se_failure = validation_sizes$se_total * (1 - unname(config$sensitivity)),
    sp_success = validation_sizes$sp_total * unname(config$specificity),
    sp_failure = validation_sizes$sp_total * (1 - unname(config$specificity))
  )
}

chapter3_config <- function(scenario) {
  registry <- chapter3_scenarios()
  if (!scenario %in% registry$scenario) {
    stop("Unknown scenario: ", scenario, ". Valid scenarios: ",
         paste(registry$scenario, collapse = ", "))
  }

  config <- canonical_config()
  if (scenario == "null") {
    config$p_y <- c(`0` = 0.35, `1` = 0.35)
  } else if (scenario == "reversal") {
    interval <- revision_induced_reversal_interval(config)
    config$p_y["1"] <- mean(interval)
  } else if (scenario == "delay_short") {
    config$central_delay_median_months <-
      0.5 * config$central_delay_median_months
    config$central_delay_q90_months <-
      0.5 * config$central_delay_q90_months
  } else if (scenario == "delay_long") {
    config$central_delay_median_months <-
      2 * config$central_delay_median_months
    config$central_delay_q90_months <-
      2 * config$central_delay_q90_months
  } else if (scenario == "transport_moderate") {
    config$sensitivity <- mix_arm_revision(config$sensitivity, 0.5)
    config$specificity <- mix_arm_revision(config$specificity, 0.5)
  } else if (scenario == "transport_severe") {
    config$sensitivity <- mix_arm_revision(config$sensitivity, 1)
    config$specificity <- mix_arm_revision(config$specificity, 1)
  } else if (scenario == "surrogate_weak") {
    config$sensitivity <- weaken_binary_accuracy(config$sensitivity)
    config$specificity <- weaken_binary_accuracy(config$specificity)
    config$historical_revision_counts <- validation_counts_at_accuracy(config)
  } else if (scenario == "surrogate_strong") {
    config$sensitivity <- strengthen_binary_accuracy(config$sensitivity)
    config$specificity <- strengthen_binary_accuracy(config$specificity)
    config$historical_revision_counts <- validation_counts_at_accuracy(config)
  }

  config$scenario <- scenario
  config
}

chapter3_scenario_parameters <- function() {
  registry <- chapter3_scenarios()
  rows <- lapply(registry$scenario, function(scenario) {
    config <- chapter3_config(scenario)
    data.frame(
      scenario = scenario,
      family = registry$family[registry$scenario == scenario],
      level = registry$level[registry$scenario == scenario],
      p0 = unname(config$p_y["0"]),
      p1 = unname(config$p_y["1"]),
      sensitivity0 = unname(config$sensitivity["0"]),
      sensitivity1 = unname(config$sensitivity["1"]),
      specificity0 = unname(config$specificity["0"]),
      specificity1 = unname(config$specificity["1"]),
      delay_median = config$central_delay_median_months,
      delay_q90 = config$central_delay_q90_months,
      n_rep = FORMAL_N_REP,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

validate_formal_chapter3_config <- function() {
  registry <- chapter3_scenarios()
  stopifnot(!anyDuplicated(registry$scenario))
  stopifnot(all(vapply(registry$scenario, function(s) {
    config <- chapter3_config(s)
    all(config$p_y > 0 & config$p_y < 1) &&
      all(config$sensitivity > 0 & config$sensitivity < 1) &&
      all(config$specificity > 0 & config$specificity < 1) &&
      config$n == config$n_months * config$patients_per_month
  }, logical(1))))
  invisible(TRUE)
}
