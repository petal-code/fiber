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
    colnames(df) <- c("weight", "R0", "prop_funeral", "p_hcw_hosp",
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
      ph_med         = round(wq(df$p_hcw_hosp,   w, 0.5), 3),
      mean_deaths    = round(sum(w * df$n_deaths)),
      mean_hcw       = round(sum(w * df$n_hcw_deaths)),
      mean_duration  = round(sum(w * df$duration)),
      mean_takeoff = round(sum(w * df$takeoff), 3)
    )
  }

  do.call(rbind, rows)
}

abc_progress(PACKAGE_ROOT, tolerance_target = 1.0)
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
