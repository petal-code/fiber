#!/usr/bin/env bash
#
# Self-contained equivalence check: proves the CURRENT working tree produces
# bit-for-bit identical() branching_process_main() output to a REFERENCE commit
# (default origin/main) across a grid of scenarios x seeds. No dev/ companion
# files needed -- the runner + comparator are embedded and used for BOTH trees,
# so only the package code differs. Uses an isolated worktree; your branch and
# working tree are never touched.
#
#   bash dev/verify_equivalence.sh [REF]      # REF defaults to origin/main
#
set -euo pipefail
command -v Rscript >/dev/null || { echo "Rscript not found on PATH"; exit 1; }
REPO="$(git rev-parse --show-toplevel)"; cd "$REPO"
REF="${1:-origin/main}"
git fetch -q origin main 2>/dev/null || echo "   (offline; using local refs)"
REF_SHA="$(git rev-parse --short "$REF")"

WORK="$(mktemp -d)"; WT="$WORK/ref"
RUNNER="$WORK/run.R"; CMP="$WORK/cmp.R"; CUR="$WORK/cur.rds"; REFR="$WORK/ref.rds"
cleanup(){ git worktree remove --force "$WT" 2>/dev/null||true; git worktree prune 2>/dev/null||true; rm -rf "$WORK"; }
trap cleanup EXIT

# ---- embedded runner (identical for both trees) ----
cat > "$RUNNER" <<'RSCRIPT'
a <- commandArgs(trailingOnly=TRUE); pkg <- a[[1]]; out <- a[[2]]
suppressMessages(devtools::load_all(pkg, quiet=TRUE))
fixed <- function(v) function(n) rep(v, n)
base_args <- list(
  mn_offspring_genPop=1.5, overdisp_offspring_genPop=0.5, Tg_shape_genPop=2, Tg_rate_genPop=0.15,
  mn_offspring_hcw=1.5, overdisp_offspring_hcw=0.5, Tg_shape_hcw=2, Tg_rate_hcw=0.15,
  mn_offspring_funeral=1.5, overdisp_offspring_funeral=0.5, Tg_shape_funeral=10, Tg_rate_funeral=5,
  incubation_period=fixed(5), onset_to_hospitalisation=fixed(3),
  onset_to_death=fixed(7), onset_to_recovery=fixed(10),
  hospitalisation_to_death=fixed(7), hospitalisation_to_recovery=fixed(10),
  prob_symptomatic=0.9, prob_hospitalised_hcw=0.5, prob_hospitalised_genPop=0.5,
  prob_death_comm=0.5, prob_death_hosp=0.3,
  prob_hcw_cond_genPop_comm=0.01, prob_hcw_cond_genPop_hospital=0.3,
  prob_hcw_cond_hcw_comm=0.01, prob_hcw_cond_hcw_hospital=0.3,
  prob_hospital_cond_hcw_preAdm=0.5, ppe_coverage_hcw=0.5, ppe_efficacy=1,
  prop_etu=1, etu_efficacy=0.5, general_hospital_quarantine_efficacy=0.5,
  p_unsafe_funeral_comm_hcw=0.5, p_unsafe_funeral_hosp_hcw=0.2,
  p_unsafe_funeral_comm_genPop=0.5, p_unsafe_funeral_hosp_genPop=0.2,
  safe_funeral_efficacy=1.0, prob_hcw_cond_funeral_hcw=0.05, prob_hcw_cond_funeral_genPop=0.05,
  population=5000, hcw_per_capita=0.02, check_final_size=1500L, seeding_cases=3, seed=1L)
scenarios <- list(
  base         = list(check_final_size=1500L),
  takeoff      = list(check_final_size=1500L, mn_offspring_genPop=3, mn_offspring_hcw=3, mn_offspring_funeral=3),
  cap_200      = list(check_final_size=200L,  mn_offspring_genPop=3, mn_offspring_hcw=3, mn_offspring_funeral=3),
  obv_on       = list(check_final_size=1500L, obv_pep_enabled=TRUE, obv_pep_coverage=0.8,
                      mn_offspring_genPop=2, mn_offspring_hcw=2, mn_offspring_funeral=2),
  time_varying = list(check_final_size=1500L,
                      mn_offspring_genPop=make_time_varying(c(0,30,80), c(1,3,0.5)),
                      mn_offspring_hcw=2, mn_offspring_funeral=2),
  susc_limited = list(check_final_size=2000L, population=1500,
                      mn_offspring_genPop=3, mn_offspring_hcw=3, mn_offspring_funeral=3))
seeds <- 1:6
res <- list()
for (sc in names(scenarios)) for (s in seeds) {
  key <- paste0(sc, "__seed", s)
  res[[key]] <- tryCatch(
    do.call(branching_process_main,
            utils::modifyList(utils::modifyList(base_args, scenarios[[sc]]), list(seed=s))),
    error=function(e) paste0("ERROR: ", conditionMessage(e)))
}
saveRDS(res, out)
cat(sprintf("   ran %d scenario x seed combinations\n", length(res)))
RSCRIPT

# ---- embedded comparator (full-object identical(); drills into tdf on mismatch) ----
cat > "$CMP" <<'RSCRIPT'
a <- commandArgs(trailingOnly=TRUE); x <- readRDS(a[[1]]); y <- readRDS(a[[2]])
keys <- union(names(x), names(y)); ndiff <- 0L
for (k in keys) {
  if (identical(x[[k]], y[[k]])) next
  ndiff <- ndiff + 1L; cat(sprintf("DIFF  %s\n", k))
  xt <- x[[k]]$tdf; yt <- y[[k]]$tdf
  if (is.data.frame(xt) && is.data.frame(yt)) {
    for (col in union(names(xt), names(yt))) {
      if (!identical(xt[[col]], yt[[col]])) {
        tag <- if (isTRUE(all.equal(xt[[col]], yt[[col]], check.attributes=FALSE)))
                 "type/attr only (values equal)" else "VALUES differ"
        cat(sprintf("      tdf$%-28s %s\n", col, tag))
      }
    }
  }
}
if (ndiff == 0L) cat(sprintf("\nPASS: all %d (scenario x seed) results identical().\n", length(keys)))
if (ndiff >  0L) { cat(sprintf("\nFAIL: %d of %d differ.\n", ndiff, length(keys))); quit(status=1) }
RSCRIPT

echo ">> current tree : $(git rev-parse --short HEAD) ($(git rev-parse --abbrev-ref HEAD))"
echo ">> reference    : $REF_SHA ($REF)"
echo ">> [1/3] running scenarios on current tree ..."
Rscript "$RUNNER" "$REPO" "$CUR"
echo ">> [2/3] running scenarios on reference (isolated worktree) ..."
git worktree add -q --detach "$WT" "$REF"
Rscript "$RUNNER" "$WT" "$REFR"
echo ">> [3/3] comparing ..."
echo "-------------------------------------------------------------"
Rscript "$CMP" "$CUR" "$REFR"
echo "-------------------------------------------------------------"
