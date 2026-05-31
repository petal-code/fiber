#!/usr/bin/env bash
## dev/benchmark_compare.sh
## ---------------------------------------------------------------------------
## One-command A/B timing: runs dev/benchmark_timing.R on a reference commit
## (the OLD model) and on the current working tree (the NEW model), back to back
## on the same machine, and prints both. Optional convenience wrapper -- you can
## equally run dev/benchmark_timing.R by hand in each checkout.
##
## Requires: git, and R with the package's dependencies + devtools installed.
##
## Usage:
##   bash dev/benchmark_compare.sh [REF]
##   REF defaults to e384a31 (the commit immediately before the optimisation).
##
## The benchmark script from the CURRENT tree is used for BOTH runs, so only the
## package code differs.
## ---------------------------------------------------------------------------
set -euo pipefail

REF="${1:-e384a318e0f09f5c0498e377db06dcd1e094ed5a}"
REPO="$(git rev-parse --show-toplevel)"
WT="$(mktemp -d)/wt"   # git worktree add needs a not-yet-existing path

cleanup() { git -C "$REPO" worktree remove --force "$WT" >/dev/null 2>&1 || true; }
trap cleanup EXIT

git -C "$REPO" worktree add --detach "$WT" "$REF" >/dev/null
cp "$REPO/dev/benchmark_timing.R" "$WT/dev/benchmark_timing.R"

echo "============================================================"
echo " OLD model  (reference: $REF)"
echo "============================================================"
( cd "$WT"   && Rscript -e 'source("dev/benchmark_timing.R")' )

echo
echo "============================================================"
echo " NEW model  (working tree: $REPO)"
echo "============================================================"
( cd "$REPO" && Rscript -e 'source("dev/benchmark_timing.R")' )
