#!/usr/bin/env bash
## dev/verify_equivalence.sh
## ---------------------------------------------------------------------------
## Proves the columnar rewrite of branching_process_main() is output-equivalent
## to a reference commit, bit-for-bit (identical()), across a grid of scenarios
## x seeds (see dev/eq_run.R). NOT part of the package.
##
## Requires: git, and R with the package's dependencies + devtools installed.
##
## Usage:
##   bash dev/verify_equivalence.sh [REF]
##
##   REF defaults to the commit immediately BEFORE the optimisation (the
##   "Add dev profiling harness" commit). Pass an explicit ref/sha to compare
##   against something else, e.g. `bash dev/verify_equivalence.sh main`.
##
## The harness scripts (eq_run.R / eq_compare.R) from the CURRENT working tree
## are used for BOTH versions, so the only thing that varies is the package code.
## ---------------------------------------------------------------------------
set -euo pipefail

REF="${1:-e384a318e0f09f5c0498e377db06dcd1e094ed5a}"
REPO="$(git rev-parse --show-toplevel)"
WT="$(mktemp -d)/wt"   # git worktree add needs a not-yet-existing path
REF_OUT="$(mktemp -u).rds"
NEW_OUT="$(mktemp -u).rds"

cleanup() {
  git -C "$REPO" worktree remove --force "$WT" >/dev/null 2>&1 || true
  rm -f "$REF_OUT" "$NEW_OUT" 2>/dev/null || true
}
trap cleanup EXIT

echo ">> reference commit : $REF"
echo ">> working tree     : $REPO"
git -C "$REPO" worktree add --detach "$WT" "$REF" >/dev/null

## Hold the harness constant across both versions.
cp "$REPO/dev/eq_run.R" "$WT/dev/eq_run.R"

echo ">> running reference ..."
( cd "$WT"   && Rscript dev/eq_run.R "$REF_OUT" )

echo ">> running working tree ..."
( cd "$REPO" && Rscript dev/eq_run.R "$NEW_OUT" )

echo ">> comparing ..."
Rscript "$REPO/dev/eq_compare.R" "$REF_OUT" "$NEW_OUT"
