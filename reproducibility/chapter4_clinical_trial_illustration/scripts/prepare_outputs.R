root <- "chapter4_clinical_trial_illustration"
methods <- c("CR", "Y-DBCD", "FL-CARA", "I-DBCD", "SP-RAR", "P-RAR", "ReRAR")
result <- readRDS(file.path(root, "results", "replicate_results.rds"))
trace <- readRDS(file.path(root, "results", "monthly_trace.rds"))
main <- read.csv(file.path(root, "tables", "table2_main_results.csv"),
                 check.names = FALSE)
if ("expected_nonresponses" %in% names(main))
  names(main)[names(main) == "expected_nonresponses"] <- "expected_progressions"
if ("nonresponses_avoided_vs_cr" %in% names(main))
  names(main)[names(main) == "nonresponses_avoided_vs_cr"] <-
    "progressions_avoided_vs_cr"

stopifnot(all(table(result$method) == 10000L))
stopifnot(setequal(unique(result$method), methods))
stopifnot(all(table(trace$method) == 150000L))

data_dir <- file.path("output", "figure_data", "clinical_trial")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
write_clean <- function(x, name) {
  write.csv(x, file.path(data_dir, name), row.names = FALSE, quote = FALSE)
}


trajectory <- read.csv(file.path(root, "results", "monthly_trajectory.csv"),
                       check.names = FALSE)
for (method in methods) {
  x <- trajectory[trajectory$method == method, ]
  write_clean(x, paste0("trajectory_", gsub("-", "_", method), ".csv"))
}


for (j in seq_along(methods)) {
  method <- methods[j]
  tag <- gsub("-", "_", method)
  x <- result$allocation_1[result$method == method]
  den <- density(x, from = .10, to = .90, n = 501, adjust = .85)
  y <- (length(methods) - j + 1) + .72 * den$y / max(den$y)
  write_clean(data.frame(x = den$x, y = y),
              paste0("allocation_ridge_", tag, ".csv"))

  est <- result$contrast_stabilized_estimate[result$method == method]
  den_est <- density(est, from = -.05, to = .45, n = 501, adjust = .90)
  y_est <- (length(methods) - j + 1) + .72 * den_est$y / max(den_est$y)
  write_clean(data.frame(x = den_est$x, y = y_est),
              paste0("estimate_ridge_", tag, ".csv"))

  benefit <- sort(result$failures_avoided_vs_cr[result$method == method])
  keep <- unique(round(seq(1, length(benefit), length.out = 601)))
  write_clean(data.frame(
    benefit = benefit[keep], probability = keep / length(benefit)
  ), paste0("benefit_ecdf_", tag, ".csv"))
}


table1 <- data.frame(
  Treatment = c("Placebo", "Cabozantinib"),
  N = c(111, 219),
  `IRC NP / Investigator NP` = c(36, 131),
  `IRC NP / Investigator PD` = c(19, 19),
  `IRC PD / Investigator NP` = c(8, 25),
  `IRC PD / Investigator PD` = c(48, 44),
  `IRC nonprogression rate` = c(55 / 111, 150 / 219),
  `Investigator nonprogression rate` = c(44 / 111, 156 / 219),
  check.names = FALSE
)
write.csv(table1, file.path(root, "tables", "table1_published_evidence.csv"),
          row.names = FALSE)
paired <- data.frame(treatment = rep(c("Placebo", "Cabozantinib"), each = 4),
                     Y = rep(c(0, 0, 1, 1), 2), Z = rep(c(0, 1), 4),
                     count = c(48, 8, 19, 36, 44, 25, 19, 131))
paired$row_percentage <- ave(paired$count, paired$treatment, paired$Y,
                             FUN = function(x) 100 * x / sum(x))
write_clean(paired, "figure6_paired_assessments.csv")
revision <- data.frame(
  treatment = table1$Treatment,
  investigator_nonprogression = table1[["Investigator nonprogression rate"]],
  irc_nonprogression = table1[["IRC nonprogression rate"]]
)
revision$net_revision_pp <- 100 * (revision$irc_nonprogression - revision$investigator_nonprogression)
revision$discordant_n <- c(27L, 44L)
revision$discordant_percent <- 100 * revision$discordant_n / table1$N
write_clean(revision, "figure6_arm_revision.csv")
write_clean(data.frame(
  assessment = c("Investigator", "IRC"),
  treatment_contrast = c(diff(revision$investigator_nonprogression), diff(revision$irc_nonprogression))
), "figure6_treatment_contrast.csv")


main <- main[match(methods, main$method), ]
digits <- c(allocation_experimental = 4, expected_progressions = 3,
            progressions_avoided_vs_cr = 3, bias = 4,
            empirical_se = 4, mean_estimated_se = 4,
            power = 4, coverage = 4)
for (nm in names(digits)) main[[nm]] <- round(main[[nm]], digits[[nm]])
write.csv(main, file.path(root, "tables", "table2_main_results.csv"),
          row.names = FALSE)


table3 <- do.call(rbind, lapply(methods, function(method) {
  x <- result[result$method == method, ]
  data.frame(
    method = method,
    allocation_q10 = unname(quantile(x$allocation_1, .10)),
    allocation_median = median(x$allocation_1),
    allocation_q90 = unname(quantile(x$allocation_1, .90)),
    benefit_q10 = unname(quantile(x$failures_avoided_vs_cr, .10)),
    benefit_median = median(x$failures_avoided_vs_cr),
    benefit_q90 = unname(quantile(x$failures_avoided_vs_cr, .90)),
    probability_benefit_positive = mean(x$failures_avoided_vs_cr > 0)
  )
}))
numeric_columns <- setdiff(names(table3), "method")
table3[numeric_columns] <- lapply(table3[numeric_columns], round, 4)
write.csv(table3, file.path(root, "tables", "table3_distribution_summary.csv"),
          row.names = FALSE)

cat("Prepared all-seven-method Chapter 4 figure data and tables.\n")
