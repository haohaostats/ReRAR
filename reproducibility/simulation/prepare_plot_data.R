
root <- normalizePath(getwd(), winslash = "/")
formal <- file.path(root, "simulation", "formal_output", "tables")
results <- file.path(root, "simulation", "formal_results")
out <- file.path(root, "output", "figure_data", "simulation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

methods <- c("CR", "Y-DBCD", "FL-CARA", "I-DBCD", "SP-RAR", "P-RAR", "ReRAR")
safe <- setNames(c("CR", "YDBCD", "FLCARA", "IDBCD", "SPRAR", "PRAR", "ReRAR"), methods)
write_clean <- function(x, name) write.table(
  x, file.path(out, name), sep = ",", row.names = FALSE,
  col.names = TRUE, quote = FALSE, na = "nan"
)

trace <- read.csv(file.path(results, "positive", "monthly_trace_summary.csv"),
                  check.names = FALSE)
for (m in methods) {
  z <- trace[trace$method == m,
             c("month", "randomization_probability_mean",
               "randomization_probability_q10", "randomization_probability_q90")]
  names(z) <- c("month", "mean", "q10", "q90")
  write_clean(z, paste0("trajectory_", safe[[m]], ".csv"))
}

positive <- read.csv(file.path(formal, "table2.csv"))
positive$method_key <- unname(safe[positive$Method])
write_clean(positive[, c("method_key", "Expected.failures", "Power")],
            "welfare_power.csv")

reversal <- read.csv(file.path(formal, "table4.csv"))
reversal$method_key <- unname(safe[reversal$Method])
cr_fail <- reversal$Expected.failures[reversal$Method == "CR"]
reversal$failures_avoided <- cr_fail - reversal$Expected.failures
write_clean(reversal[, c("method_key", "Allocation.to.truly.better.arm",
                         "failures_avoided")], "reversal.csv")

sens <- read.csv(file.path(formal, "table5.csv"))
make_sensitivity <- function(name, scenarios, labels) {
  ans <- data.frame(x = seq_along(scenarios), label = labels)
  for (m in methods) {
    rows <- lapply(scenarios, function(s) {
      sens[sens$scenario == s & sens$method == m, , drop = FALSE]
    })
    ans[[safe[[m]]]] <- vapply(rows, function(z) z$expected_failures, numeric(1))
    ans[[paste0(safe[[m]], "Power")]] <- vapply(rows, function(z) z$power, numeric(1))
  }
  write_clean(ans, paste0("sensitivity_", name, ".csv"))
}

make_sensitivity("delay", c("delay_short", "positive", "delay_long"),
                 c("Short", "Reference", "Long"))
make_sensitivity("transport", c("positive", "transport_moderate", "transport_severe"),
                 c("Correct", "Moderate", "Severe"))
make_sensitivity("surrogate", c("surrogate_weak", "positive", "surrogate_strong"),
                 c("Weak", "Reference", "Strong"))


interval_rows <- function(values, y) {
  q <- unname(quantile(values, c(.10, .50, .90), na.rm = TRUE))
  data.frame(y = y, median = q[2], errminus = q[2] - q[1],
             errplus = q[3] - q[2])
}
for (scenario in c("positive", "null", "reversal")) {
  x <- readRDS(file.path(results, scenario, "replications.rds"))
  rows <- do.call(rbind, lapply(seq_along(methods), function(i) {
    m <- methods[i]
    cbind(method = safe[[m]],
          interval_rows(x$allocation_1[x$method == m], length(methods) + 1L - i))
  }))
  write_clean(rows[rows$method != "ReRAR", ],
              paste0("allocation_", scenario, "_comparators.csv"))
  write_clean(rows[rows$method == "ReRAR", ],
              paste0("allocation_", scenario, "_rerar.csv"))
}


for (scenario in c("positive", "reversal")) {
  x <- readRDS(file.path(results, scenario, "replications.rds"))
  cr <- x[x$method == "CR", c("replicate", "expected_failures")]
  names(cr)[2] <- "cr"
  adaptive <- methods[-1]
  rows <- do.call(rbind, lapply(seq_along(adaptive), function(i) {
    m <- adaptive[i]
    z <- merge(x[x$method == m, c("replicate", "expected_failures")],
               cr, by = "replicate")
    cbind(method = safe[[m]],
          interval_rows(z$cr - z$expected_failures, length(adaptive) + 1L - i))
  }))
  write_clean(rows[rows$method != "ReRAR", ],
              paste0("paired_", scenario, "_comparators.csv"))
  write_clean(rows[rows$method == "ReRAR", ],
              paste0("paired_", scenario, "_rerar.csv"))


  for (m in adaptive) {
    z <- merge(x[x$method == m, c("replicate", "expected_failures")],
               cr, by = "replicate")
    benefit <- sort(z$cr - z$expected_failures)
    ecdf_full <- data.frame(
      benefit = benefit,
      cumulative_probability = seq_along(benefit) / length(benefit)
    )
    write_clean(ecdf_full,
                paste0("paired_ecdf_", scenario, "_", safe[[m]], ".csv"))



    keep <- unique(round(seq(1, length(benefit), length.out = 401)))
    ecdf_grid <- data.frame(
      benefit = benefit[keep],
      cumulative_probability = keep / length(benefit)
    )
    write_clean(ecdf_grid,
                paste0("paired_ecdf_", scenario, "_", safe[[m]], "_plot.csv"))
  }
}


inference <- read.csv(file.path(formal, "table7.csv"))
for (scenario in c("positive", "null", "reversal")) {
  z <- inference[inference$scenario == scenario, ]
  z <- z[match(methods, z$method), ]
  z$x <- seq_along(methods)
  keep <- c("x", "method", "bias", "empirical_se", "mean_estimated_se")
  write_clean(z[z$method != "ReRAR", keep],
              paste0("inference_", scenario, "_comparators.csv"))
  write_clean(z[z$method == "ReRAR", keep],
              paste0("inference_", scenario, "_rerar.csv"))
}



positive_rep <- readRDS(file.path(results, "positive", "replications.rds"))
positive_cr <- positive_rep[positive_rep$method == "CR",
                            c("replicate", "expected_failures")]
names(positive_cr)[2] <- "cr"
thresholds <- seq(-8, 10, by = .25)
joint_curve <- data.frame(threshold = thresholds)
for (m in methods) {
  z <- merge(positive_rep[positive_rep$method == m,
                          c("replicate", "expected_failures",
                            "contrast_stabilized_reject")],
             positive_cr, by = "replicate")
  benefit <- z$cr - z$expected_failures
  joint_curve[[safe[[m]]]] <- vapply(
    thresholds,
    function(d) mean(benefit >= d & z$contrast_stabilized_reject == 1L),
    numeric(1)
  )
}
write_clean(joint_curve, "joint_welfare_efficacy.csv")


reversal_rep <- readRDS(file.path(results, "reversal", "replications.rds"))
reversal_cr <- reversal_rep[reversal_rep$method == "CR",
                            c("replicate", "expected_failures")]
names(reversal_cr)[2] <- "cr"
ellipse_points <- function(x, y, probability, n = 121L) {
  centre <- c(mean(x), mean(y))
  sigma <- stats::cov(cbind(x, y))
  root <- chol(sigma + diag(c(1e-10, 1e-10)))
  angle <- seq(0, 2 * pi, length.out = n)
  unit <- rbind(cos(angle), sin(angle))
  points <- t(centre + sqrt(stats::qchisq(probability, df = 2)) *
                t(root) %*% unit)
  data.frame(x = points[, 1], y = points[, 2])
}
reversal_means <- list()
for (m in methods) {
  z <- merge(reversal_rep[reversal_rep$method == m,
                          c("replicate", "allocation_1", "expected_failures")],
             reversal_cr, by = "replicate")
  better <- 1 - z$allocation_1
  benefit <- z$cr - z$expected_failures
  reversal_means[[m]] <- data.frame(
    method = safe[[m]], x = mean(better), y = mean(benefit)
  )
  if (m != "CR") {
    write_clean(ellipse_points(better, benefit, .50),
                paste0("reversal_ellipse_", safe[[m]], "_50.csv"))
    write_clean(ellipse_points(better, benefit, .90),
                paste0("reversal_ellipse_", safe[[m]], "_90.csv"))
  }
}
write_clean(do.call(rbind, reversal_means), "reversal_means.csv")


principal <- list(positive = positive_rep,
                  null = readRDS(file.path(results, "null", "replications.rds")),
                  reversal = reversal_rep)
for (scenario in names(principal)) {
  x <- principal[[scenario]]
  for (i in seq_along(methods)) {
    m <- methods[i]
    value <- if (scenario == "reversal") {
      1 - x$allocation_1[x$method == m]
    } else x$allocation_1[x$method == m]
    den <- density(value, from = .30, to = .70, n = 301, bw = "nrd0")
    baseline <- length(methods) + 1L - i
    ridge <- data.frame(
      x = den$x, baseline = baseline,
      curve = baseline + .72 * den$y / max(den$y)
    )
    write_clean(ridge, paste0("density_", scenario, "_", safe[[m]], ".csv"))
  }
}




pairwise_probability <- function(x, row_method, column_method) {
  a <- x[x$method == row_method, c("replicate", "expected_failures")]
  b <- x[x$method == column_method, c("replicate", "expected_failures")]
  names(a)[2] <- "a"; names(b)[2] <- "b"
  z <- merge(a, b, by = "replicate")
  mean(z$a < z$b) + .5 * mean(z$a == z$b)
}
matrix_rows <- do.call(rbind, lapply(seq_along(methods), function(i) {
  do.call(rbind, lapply(seq_along(methods), function(j) {
    if (i == j) {
      probability <- .5
    } else if (i > j) {
      probability <- pairwise_probability(positive_rep, methods[i], methods[j])
    } else {
      probability <- pairwise_probability(reversal_rep, methods[i], methods[j])
    }
    data.frame(x = j, y = length(methods) + 1L - i,
               probability = probability,
               label = sprintf("%.2f", probability))
  }))
}))
write_clean(matrix_rows, "pairwise_dominance_matrix.csv")

message("Figure data written to ", out)
