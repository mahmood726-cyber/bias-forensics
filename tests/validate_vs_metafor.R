# R cross-validation: compare Python Egger/Begg/TrimFill against metafor
# Run: Rscript tests/validate_vs_metafor.R

library(metafor)

py_results <- read.csv("C:/BiasForensics/data/output/bias_forensics_results.csv",
                        stringsAsFactors = FALSE, fileEncoding = "UTF-8")

pairwise_dir <- "C:/Models/Pairwise70/data"

# Sample 10 reviews
set.seed(42)
sample_ids <- py_results$review_id[sample(nrow(py_results), min(10, nrow(py_results)))]

cat("=== R Cross-Validation: Python vs metafor ===\n\n")

for (rid in sample_ids) {
  rda_files <- list.files(pairwise_dir, pattern = paste0("^", rid, "_"), full.names = TRUE)
  if (length(rda_files) == 0) next

  env <- new.env()
  load(rda_files[1], envir = env)
  df <- get(ls(env)[1], envir = env)

  # Select primary analysis (largest k)
  analyses <- aggregate(Study ~ Analysis.group + Analysis.number, data = df, FUN = length)
  names(analyses)[3] <- "k"
  best <- analyses[which.max(analyses$k), ]
  primary <- df[df$Analysis.group == best$Analysis.group &
                df$Analysis.number == best$Analysis.number, ]

  primary <- primary[!is.na(primary$Mean) & primary$Mean > 0 &
                     !is.na(primary$CI.start) & primary$CI.start > 0 &
                     !is.na(primary$CI.end) & primary$CI.end > 0, ]
  if (nrow(primary) < 5) next

  yi <- log(primary$Mean)
  sei <- (log(primary$CI.end) - log(primary$CI.start)) / (2 * 1.96)
  valid <- !is.na(yi) & !is.na(sei) & sei > 0
  yi <- yi[valid]; sei <- sei[valid]
  if (length(yi) < 5) next

  py_row <- py_results[py_results$review_id == rid, ]
  if (nrow(py_row) == 0) next

  tryCatch({
    fit <- rma(yi = yi, sei = sei, method = "DL")

    # Egger
    reg <- regtest(fit, model = "lm")
    r_egger_p <- reg$pval

    # Begg
    rank <- ranktest(fit)
    r_begg_p <- rank$pval

    # Trim-fill
    tf <- trimfill(fit)
    r_tf_k0 <- tf$k0
    r_tf_theta <- as.numeric(coef(tf))

    cat(sprintf("  %s (k=%d):\n", rid, length(yi)))
    cat(sprintf("    Egger:  R p=%.4f, Py p=%.4f, diff=%.4f\n",
                r_egger_p, as.numeric(py_row$egger_p), abs(r_egger_p - as.numeric(py_row$egger_p))))
    cat(sprintf("    Begg:   R p=%.4f, Py p=%.4f, diff=%.4f\n",
                r_begg_p, as.numeric(py_row$begg_p), abs(r_begg_p - as.numeric(py_row$begg_p))))
    cat(sprintf("    TF k0:  R=%d, Py=%d\n", r_tf_k0, as.integer(py_row$tf_k0)))
    cat(sprintf("    TF adj: R=%.4f, Py=%.4f, diff=%.4f\n",
                r_tf_theta, as.numeric(py_row$tf_theta_adj),
                abs(r_tf_theta - as.numeric(py_row$tf_theta_adj))))
    cat("\n")
  }, error = function(e) {
    cat(sprintf("  %s: ERROR - %s\n\n", rid, conditionMessage(e)))
  })
}
