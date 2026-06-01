## dev/eq_compare.R
## ---------------------------------------------------------------------------
## Equivalence-harness COMPARATOR. NOT part of the package.
##
## Compares two .rds files produced by dev/eq_run.R and reports, per
## (scenario x seed), whether the full returned object is identical(). On any
## mismatch it drills into the $tdf data frame column-by-column and says whether
## the difference is in VALUES or only in storage TYPE/attributes -- the latter
## being scientifically irrelevant but flagged so it can be made bit-identical.
##
## Usage: Rscript dev/eq_compare.R ref.rds new.rds
## Exit status 0 == all identical, 1 == at least one difference.
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: Rscript dev/eq_compare.R ref.rds new.rds")
a <- readRDS(args[[1L]])
b <- readRDS(args[[2L]])

keys <- union(names(a), names(b))
all_ok <- TRUE

describe_df_diff <- function(da, db) {
  if (!is.data.frame(da) || !is.data.frame(db)) {
    cat("    (one side is not a data.frame: ", class(da)[1], " vs ",
        class(db)[1], ")\n", sep = "")
    return(invisible())
  }
  if (!identical(dim(da), dim(db))) {
    cat(sprintf("    dim differs: %s vs %s\n",
                paste(dim(da), collapse = "x"), paste(dim(db), collapse = "x")))
  }
  cols <- union(names(da), names(db))
  for (col in cols) {
    ca <- da[[col]]; cb <- db[[col]]
    if (identical(ca, cb)) next
    tag <- if (is.null(ca) || is.null(cb)) {
      "MISSING COLUMN"
    } else if (!isTRUE(all.equal(ca, cb, check.attributes = FALSE))) {
      "VALUES differ"
    } else if (!identical(typeof(ca), typeof(cb))) {
      sprintf("type only: %s vs %s (values equal)", typeof(ca), typeof(cb))
    } else {
      "attributes only (values equal)"
    }
    cat(sprintf("    col %-32s %s\n", col, tag))
  }
}

for (k in keys) {
  if (identical(a[[k]], b[[k]])) next
  all_ok <- FALSE
  cat(sprintf("DIFF  %s\n", k))
  ea <- a[[k]]; eb <- b[[k]]
  if (is.character(ea) || is.character(eb)) {
    cat("    ref: ", if (is.character(ea)) ea else "<object>", "\n")
    cat("    new: ", if (is.character(eb)) eb else "<object>", "\n")
    next
  }
  # compare $tdf and the rest of the list separately
  describe_df_diff(ea$tdf, eb$tdf)
  rest_a <- ea[setdiff(names(ea), "tdf")]
  rest_b <- eb[setdiff(names(eb), "tdf")]
  if (!identical(rest_a, rest_b)) cat("    sim_info / attrs differ\n")
  # attributes on tdf (hcw_total, obv_pep_num_treated, ...)
  at_a <- attributes(ea$tdf); at_b <- attributes(eb$tdf)
  for (an in union(names(at_a), names(at_b))) {
    if (an %in% c("names", "row.names", "class")) next
    if (!identical(at_a[[an]], at_b[[an]])) cat(sprintf("    attr(tdf,'%s') differs\n", an))
  }
}

if (all_ok) {
  cat(sprintf("\nPASS: all %d (scenario x seed) results are identical().\n", length(keys)))
  quit(status = 0)
} else {
  cat("\nFAIL: differences reported above.\n")
  quit(status = 1)
}
