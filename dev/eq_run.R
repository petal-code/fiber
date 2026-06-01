## dev/eq_run.R
## ---------------------------------------------------------------------------
## Equivalence-harness RUNNER. NOT part of the package.
##
## Runs branching_process_main() over a grid of (scenario x seed) and saves the
## FULL returned objects (tdf + all attributes + sim_info) to an .rds so two
## package versions can be compared bit-for-bit with identical().
##
## This file is held CONSTANT across the two versions being compared
## (verify_equivalence.sh copies it into the reference worktree), so any
## difference in the saved output is attributable to the package code alone.
##
## Usage (normally driven by dev/verify_equivalence.sh):
##   Rscript dev/eq_run.R /path/to/out.rds
## ---------------------------------------------------------------------------

args     <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1L) args[[1L]] else "/tmp/eq_out.rds"

suppressMessages(devtools::load_all(".", quiet = TRUE))
source("dev/profile_run.R")   # provides run_sim() with a full, valid parameter set

## Scenarios = overrides passed to run_sim(). Chosen to exercise every code path
## the columnar rewrite touches:
##   * fizzle (base)            -> small N, many 0-offspring parents (else branch,
##                                 the n_offspring integer->double promotion)
##   * takeoff / caps           -> the append at scale, the n_filled counter, the
##                                 `n_filled <= check_final_size` cap guard, and
##                                 over-cap vector extension (cap200)
##   * obv_on                   -> populates obv_pep_* columns in appended rows
##   * time_varying             -> resolve-at-time transmissibility path
##   * susc_limited             -> loop exit via `susc > 0` rather than the cap
##
## Caps are kept small (<=2000) deliberately: equivalence is scale-invariant
## (the append / n_filled counter / over-cap extension logic is identical at any
## size), so small caps exercise every path while keeping the whole harness to
## ~1-2 minutes instead of tens of minutes.
scenarios <- list(
  base         = list(check_final_size = 1500L),
  takeoff_1500 = list(check_final_size = 1500L, mn_offspring_genPop = 3,
                      mn_offspring_hcw = 3, mn_offspring_funeral = 3),
  cap_200      = list(check_final_size = 200L,  mn_offspring_genPop = 3,
                      mn_offspring_hcw = 3, mn_offspring_funeral = 3),
  obv_on       = list(check_final_size = 1500L, obv_pep_enabled = TRUE,
                      obv_pep_coverage = 0.8, mn_offspring_genPop = 2,
                      mn_offspring_hcw = 2, mn_offspring_funeral = 2),
  time_varying = list(check_final_size = 1500L,
                      mn_offspring_genPop = make_time_varying(c(0, 30, 80),
                                                              c(1, 3, 0.5)),
                      mn_offspring_hcw = 2, mn_offspring_funeral = 2),
  susc_limited = list(check_final_size = 2000L, population = 1500,
                      mn_offspring_genPop = 3, mn_offspring_hcw = 3,
                      mn_offspring_funeral = 3)
)
seeds <- 1:6

out <- list()
for (sc in names(scenarios)) {
  for (s in seeds) {
    key <- paste0(sc, "__seed", s)
    out[[key]] <- tryCatch(
      do.call(run_sim, c(list(seed = s), scenarios[[sc]])),
      error = function(e) paste0("ERROR: ", conditionMessage(e))
    )
  }
}

saveRDS(out, out_path)
cat(sprintf("eq_run.R: wrote %d results (%d scenarios x %d seeds) to %s\n",
            length(out), length(scenarios), length(seeds), out_path))
