PACKAGE_ROOT <- "C:/Users/PETAL_WS_2/Documents/fiber"   # or whichever
list.files(PACKAGE_ROOT, pattern = "^(output_step|n_simul_tot_step|tolerance_step)")

step1 <- read.table(file.path(PACKAGE_ROOT, "output_step1"), header = FALSE)
colnames(step1) <- c("weight", "R0", "prop_funeral", "p_hcw_hosp",
                     "takeoff", "n_deaths", "n_hcw_deaths", "duration")
head(step1)
nrow(step1)                                   # should be 100

# How many sims did step 1 cost?
readLines(file.path(PACKAGE_ROOT, "n_simul_tot_step1"))

# Quick read on whether the particles are doing anything sensible
summary(step1)
plot(step1$R0, step1$n_deaths, log = "y",
     xlab = "R0", ylab = "mean n_deaths (took-off reps)")
abline(h = 11325, col = "red", lty = 2)   # observed

summary(step1)
plot(step1$R0, step1$n_hcw_deaths, log = "y",
     xlab = "R0", ylab = "mean n_deaths (took-off reps)")
abline(h = 513, col = "red", lty = 2)   # observed

abc_progress <- function(dir = getwd(), tolerance_target = 1.0) {

  # Helper: extract step number from filename for proper numeric sorting
  step_num <- function(f) as.integer(sub(".*_step([0-9]+)$", "\\1", f))
  sort_by_step <- function(f) f[order(step_num(f))]

  out_files <- sort_by_step(list.files(dir, pattern = "^output_step[0-9]+$",       full.names = TRUE))
  tol_files <- sort_by_step(list.files(dir, pattern = "^tolerance_step[0-9]+$",    full.names = TRUE))
  sim_files <- sort_by_step(list.files(dir, pattern = "^n_simul_tot_step[0-9]+$",  full.names = TRUE))

  if (length(out_files) == 0L) {
    cat("No steps completed yet.\n")
    return(invisible(NULL))
  }

  cat(sprintf("Steps completed: %d\n", length(out_files)))

  # Cumulative sim count from the latest file
  if (length(sim_files) > 0L) {
    cum_sims <- as.numeric(readLines(sim_files[length(sim_files)]))
    cat(sprintf("Total simulations so far: %s\n", format(cum_sims, big.mark = ",")))
  }

  # Tolerance trajectory (only present from step 2 onward)
  if (length(tol_files) > 0L) {
    tols <- vapply(tol_files, function(f) as.numeric(readLines(f)), numeric(1))
    cat("\nTolerance trajectory:\n")
    for (i in seq_along(tols)) {
      cat(sprintf("  After step %d: %.3f\n", step_num(tol_files)[i], tols[i]))
    }
    cat(sprintf("\nTarget tolerance: %.3f\n", tolerance_target))

    current <- tols[length(tols)]
    if (current <= tolerance_target) {
      cat("Reached target — algorithm should terminate after this step.\n")
    } else if (length(tols) >= 2L) {
      # Project remaining steps from recent shrinkage rate
      ratio <- tols[length(tols)] / tols[length(tols) - 1L]
      if (ratio < 1) {
        remaining <- ceiling(log(tolerance_target / current) / log(ratio))
        cat(sprintf("Projected remaining steps: ~%d (assuming current shrinkage rate continues)\n",
                    remaining))
      } else {
        cat("Tolerance is not decreasing — the algorithm may be stuck.\n")
      }
    }
  }

  # Particle cloud from the latest step
  last_out <- read.table(out_files[length(out_files)], header = FALSE)
  colnames(last_out) <- c("weight", "R0", "prop_funeral", "p_hcw_hosp",
                          "takeoff", "n_deaths", "n_hcw_deaths", "duration")
  cat(sprintf("\nLatest particle cloud (step %d):\n", step_num(out_files)[length(out_files)]))
  cat(sprintf("  R0:           median = %.3f, 95%% CI = [%.3f, %.3f]\n",
              median(last_out$R0),
              quantile(last_out$R0, 0.025),
              quantile(last_out$R0, 0.975)))
  cat(sprintf("  prop_funeral: median = %.3f, 95%% CI = [%.3f, %.3f]\n",
              median(last_out$prop_funeral),
              quantile(last_out$prop_funeral, 0.025),
              quantile(last_out$prop_funeral, 0.975)))
  cat(sprintf("  p_hcw_hosp:   median = %.3f, 95%% CI = [%.3f, %.3f]\n",
              median(last_out$p_hcw_hosp),
              quantile(last_out$p_hcw_hosp, 0.025),
              quantile(last_out$p_hcw_hosp, 0.975)))

  invisible(last_out)
}

abc_compare_steps <- function(dir = getwd()) {

  step_num <- function(f) as.integer(sub(".*_step([0-9]+)$", "\\1", f))
  sort_by_step <- function(f) f[order(step_num(f))]

  out_files <- sort_by_step(list.files(dir, pattern = "^output_step[0-9]+$",      full.names = TRUE))
  tol_files <- sort_by_step(list.files(dir, pattern = "^tolerance_step[0-9]+$",   full.names = TRUE))
  sim_files <- sort_by_step(list.files(dir, pattern = "^n_simul_tot_step[0-9]+$", full.names = TRUE))

  if (length(out_files) == 0L) { cat("No steps completed yet.\n"); return(invisible(NULL)) }

  # Weighted quantile helper
  wq <- function(x, w, p) {
    ord <- order(x); xs <- x[ord]; ws <- w[ord] / sum(w)
    approx(cumsum(ws), xs, p, rule = 2, ties = "ordered")$y
  }

  tol_lookup <- if (length(tol_files) > 0L)
    setNames(vapply(tol_files, function(f) as.numeric(readLines(f)), numeric(1)),
             as.character(step_num(tol_files))) else c()
  cum_lookup <- if (length(sim_files) > 0L)
    setNames(vapply(sim_files, function(f) as.numeric(readLines(f)), numeric(1)),
             as.character(step_num(sim_files))) else c()

  rows <- list(); prev_cum <- 0
  for (f in out_files) {
    s <- step_num(f)
    df <- read.table(f, header = FALSE)
    colnames(df) <- c("weight", "R0", "prop_funeral", "hcw_risk_scalar",
                      "takeoff", "n_deaths", "n_hcw_deaths", "duration")
    w <- df$weight / sum(df$weight)

    cum_s     <- if (as.character(s) %in% names(cum_lookup)) cum_lookup[[as.character(s)]] else NA_real_
    sims_step <- if (!is.na(cum_s)) cum_s - prev_cum else NA_real_
    if (!is.na(cum_s)) prev_cum <- cum_s

    rows[[length(rows) + 1L]] <- data.frame(
      step           = s,
      tol            = if (as.character(s) %in% names(tol_lookup)) round(tol_lookup[[as.character(s)]], 3) else NA_real_,
      sims_this_step = sims_step,
      cum_sims       = cum_s,
      ESS            = round(1 / sum(w^2), 1),
      R0_med         = round(wq(df$R0,           w, 0.5), 3),
      R0_lo          = round(wq(df$R0,           w, 0.025), 3),
      R0_hi          = round(wq(df$R0,           w, 0.975), 3),
      pf_med         = round(wq(df$prop_funeral, w, 0.5), 3),
      h_scalar_med   = round(wq(df$p_hcw_hosp,   w, 0.5), 3),
      mean_deaths    = round(sum(w * df$n_deaths)),
      mean_hcw       = round(sum(w * df$n_hcw_deaths)),
      mean_duration  = round(sum(w * df$duration)),
      mean_takeoff = round(sum(w * df$takeoff), 3)
    )
  }

  do.call(rbind, rows)
}

abc_progress(PACKAGE_ROOT, tolerance_target = 0.5)
abc_compare_steps(PACKAGE_ROOT)

step3 <- read.table(file.path(PACKAGE_ROOT, "output_step3"), header = FALSE)
colnames(step3) <- c("weight", "R0", "prop_funeral", "p_hcw_hosp",
                     "takeoff", "n_deaths", "n_hcw_deaths", "duration")
summary(step3$takeoff)
table(step3$takeoff)   # cleaner — counts of distinct values

step_files <- list.files(PACKAGE_ROOT, pattern = "^output_step[0-9]+$", full.names = TRUE)
step_times <- file.info(step_files)$mtime
data.frame(file = basename(step_files),
           finished = step_times,
           gap = c(NA, diff(step_times)))

reconstruct_abc_result <- function(dir = getwd(), step = NULL) {

  step_num <- function(f) as.integer(sub(".*_step([0-9]+)$", "\\1", f))

  out_files <- list.files(dir, pattern = "^output_step[0-9]+$", full.names = TRUE)
  out_files <- out_files[order(step_num(out_files))]

  if (length(out_files) == 0L) stop("No output_step files found in ", dir)

  # Use the latest completed step unless user specifies one
  if (is.null(step)) {
    target_file <- out_files[length(out_files)]
    step <- step_num(target_file)
  } else {
    target_file <- file.path(dir, paste0("output_step", step))
    if (!file.exists(target_file)) stop("output_step", step, " not found")
  }

  message("Reconstructing result from step ", step)

  # Read the particle population
  df <- read.table(target_file, header = FALSE)
  colnames(df) <- c("weight", "R0", "prop_funeral", "p_hcw_hosp",
                    "takeoff", "n_deaths", "n_hcw_deaths", "duration")

  # Get tolerance for this step (NA if step 1)
  tol_file <- file.path(dir, paste0("tolerance_step", step))
  epsilon <- if (file.exists(tol_file)) as.numeric(readLines(tol_file)) else NA_real_

  # Get cumulative simulation count
  sim_file <- file.path(dir, paste0("n_simul_tot_step", step))
  nsim <- if (file.exists(sim_file)) as.numeric(readLines(sim_file)) else NA_real_

  # Recompute normalization SDs from step 1 (matches what EasyABC uses)
  step1_file <- file.path(dir, "output_step1")
  if (file.exists(step1_file)) {
    step1 <- read.table(step1_file, header = FALSE)
    stats_norm <- apply(step1[, 5:8], 2, sd)
  } else {
    stats_norm <- apply(df[, 5:8], 2, sd)
  }

  list(
    param               = as.matrix(df[, c("R0", "prop_funeral", "p_hcw_hosp")]),
    stats               = as.matrix(df[, c("takeoff", "n_deaths", "n_hcw_deaths", "duration")]),
    weights             = df$weight / sum(df$weight),
    stats_normalization = as.numeric(stats_norm),
    epsilon             = epsilon,
    nsim                = nsim,
    computime           = NA_real_,           # not recoverable from files
    step_reconstructed_from = step
  )
}

# Reconstruct from the latest completed step (step 5 in your case)
result <- reconstruct_abc_result(PACKAGE_ROOT)

# Or specify a particular step
result <- reconstruct_abc_result(PACKAGE_ROOT, step = 4)

# Now use it the same way as a real ABC_sequential result
posterior <- as.data.frame(result$param)
colnames(posterior) <- c("R0", "prop_funeral", "prob_hcw_hosp")

print(apply(posterior, 2, quantile, probs = c(0.025, 0.5, 0.975)))

# All the downstream plotting code from section 10 works unchanged
par(mfrow = c(1, 3))
for (j in seq_len(ncol(posterior))) {
  hist(posterior[, j], breaks = 30, main = colnames(posterior)[j])
  abline(v = quantile(posterior[, j], c(0.025, 0.5, 0.975)),
         lty = c(2, 1, 2), col = "red")
}
par(mfrow = c(1, 1))
pairs(posterior, pch = 16, cex = 0.5, col = adjustcolor("steelblue", alpha = 0.4))



# === Posterior-predictive check ===

sim_stats <- as.data.frame(result$stats)
colnames(sim_stats) <- c("takeoff", "n_deaths", "n_hcw_deaths", "duration")

# Resample by weights so the histogram reflects the weighted posterior
set.seed(1)
idx <- sample(seq_len(nrow(sim_stats)),
              size    = 10000,
              replace = TRUE,
              prob    = result$weights)
sim_stats_post <- sim_stats[idx, ]

# Observed values (matches your script)
observed <- c(takeoff      = 1.0,
              n_deaths     = 11325,
              n_hcw_deaths = 513,
              duration     = 365)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (s in names(observed)) {
  x  <- sim_stats_post[[s]]
  qs <- quantile(x, probs = c(0.025, 0.5, 0.975))

  hist(x,
       breaks = 10,
       main   = paste0("Posterior predictive: ", s),
       xlab   = s,
       col    = adjustcolor("steelblue", alpha = 0.6),
       border = "white")

  # Posterior median (solid) and 95% CI (dashed)
  abline(v = qs, col = "darkblue", lty = c(2, 1, 2), lwd = c(1, 2, 1))

  # Observed value
  abline(v = observed[s], col = "red", lwd = 2.5)

  legend("topleft",
         legend = c(paste0("Observed: ",  signif(observed[s], 3)),
                    paste0("Median: ",    signif(qs[2], 3)),
                    paste0("95% CI: [",   signif(qs[1], 3), ", ", signif(qs[3], 3), "]")),
         col    = c("red", "darkblue", NA),
         lty    = c(1, 1, NA),
         lwd    = c(2.5, 2, NA),
         bty    = "n",
         cex    = 0.75)
}
par(mfrow = c(1, 1))

# Single-panel summary: simulated vs observed, log scale where useful
par(mfrow = c(1, 1))
plot(NA, xlim = c(0.5, 4.5), ylim = c(0, 1.5),
     xaxt = "n", xlab = "", ylab = "Simulated / Observed",
     main = "Posterior-predictive fit ratio")
axis(1, at = 1:4, labels = names(observed))
abline(h = 1, lty = 2, col = "red")
for (i in seq_along(observed)) {
  x   <- sim_stats_post[[names(observed)[i]]] / observed[names(observed)[i]]
  qs  <- quantile(x, c(0.025, 0.5, 0.975))
  segments(i, qs[1], i, qs[3], lwd = 2, col = "darkblue")
  points(i, qs[2], pch = 16, cex = 1.5, col = "darkblue")
}
