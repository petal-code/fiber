## Contact-first transmission: risk-tier structure and R0 <-> baseline-risk algebra.
##
## The model can generate transmission in one of two ways (see the offspring
## functions and branching_process_main):
##
##   * "offspring mode" (legacy): the Negative Binomial draw IS the number of
##     secondary infections. Supply `mn_offspring_*` / `overdisp_offspring_*`.
##   * "contact mode": the Negative Binomial draw is the number of CONTACTS
##     (exposure events). Each contact is assigned a risk tier, and transmits
##     with probability `baseline_risk * relative_risk[tier]`. Supply
##     `mn_contacts_*` / `overdisp_contacts_*` plus either `baseline_risk_*`
##     or a target `r0_*` to solve the baseline risk from.
##
## Because each contact's tier is drawn independently, a contact's marginal
## transmission probability is `baseline_risk * mean_relative_risk` regardless
## of which tier it lands in. That makes the expected offspring count
##
##     R0 = mn_contacts * baseline_risk * mean_relative_risk
##
## which is the identity that all of the helpers below invert, and it is why
## the contact overdispersion never enters R0: thinning a Negative Binomial
## does not change its mean.

#' Define a contact risk-tier structure
#'
#' Builds the risk-tier structure used by fiber's contact-first transmission
#' model. Contacts are partitioned into tiers by \code{fractions}; each tier
#' carries a \code{relative_risk} multiplier applied to a single baseline
#' per-contact transmission probability. The tier named by \code{reference} is
#' the one the baseline risk refers to: all relative risks are divided through
#' by that tier's value, so the baseline risk parameter is literally the
#' per-contact transmission probability of the reference tier.
#'
#' The default is five equally sized tiers with equal relative risk, i.e. a
#' deliberately flat structure that adds the contact bookkeeping without
#' imposing any heterogeneity. Use \code{\link{contact_risk_gradient}} to build
#' a graded structure from a single spread parameter.
#'
#' @param fractions Numeric vector of tier shares. Must be non-negative and sum
#'   to 1. Its length sets the number of tiers. Defaults to five equal shares.
#' @param relative_risk Numeric vector, same length as \code{fractions}. Risk of
#'   each tier relative to the reference tier. Must be positive and finite.
#' @param reference Integer index (1-based) of the reference tier -- the tier the
#'   baseline risk parameter describes. Defaults to the first tier.
#' @param labels Optional character vector of tier names, same length as
#'   \code{fractions}. Used in the contact log and summaries. Defaults to
#'   \code{"risk_1"}, \code{"risk_2"}, ...
#'
#' @return An object of class \code{"fiber_contact_risk"}: a list with the
#'   validated \code{fractions}, normalised \code{relative_risk}, \code{reference},
#'   \code{labels}, \code{n_levels}, and the derived \code{mean_relative_risk}
#'   (the fraction-weighted mean, which scales R0) and \code{max_relative_risk}
#'   (which caps the feasible baseline risk at \code{1 / max_relative_risk}).
#'
#' @examples
#' ## Default: five flat tiers
#' make_contact_risk()
#'
#' ## Most contacts low risk, a small high-risk tail
#' make_contact_risk(fractions     = c(0.6, 0.25, 0.1, 0.04, 0.01),
#'                   relative_risk = c(1, 2, 5, 10, 25))
#' @export
make_contact_risk <- function(fractions     = rep(0.2, 5),
                              relative_risk = rep(1, 5),
                              reference     = 1L,
                              labels        = NULL) {

  if (!is.numeric(fractions) || length(fractions) < 1L ||
      any(!is.finite(fractions)) || any(fractions < 0)) {
    stop("`fractions` must be a non-empty numeric vector of finite, non-negative values.",
         call. = FALSE)
  }
  if (abs(sum(fractions) - 1) > 1e-8) {
    stop(sprintf("`fractions` must sum to 1 (got %s).", format(sum(fractions))),
         call. = FALSE)
  }

  n_levels <- length(fractions)

  if (!is.numeric(relative_risk) || length(relative_risk) != n_levels ||
      any(!is.finite(relative_risk)) || any(relative_risk <= 0)) {
    stop(sprintf("`relative_risk` must be a numeric vector of %d finite, strictly positive values (same length as `fractions`).",
                 n_levels), call. = FALSE)
  }

  if (!is.numeric(reference) || length(reference) != 1L || is.na(reference) ||
      reference != round(reference) || reference < 1 || reference > n_levels) {
    stop(sprintf("`reference` must be a single integer index between 1 and %d.", n_levels),
         call. = FALSE)
  }
  reference <- as.integer(reference)

  ## Normalise so the reference tier has relative risk exactly 1. The baseline
  ## risk parameter is then that tier's own per-contact transmission probability,
  ## whether or not the user happened to write its relative risk as 1.
  relative_risk <- relative_risk / relative_risk[reference]

  if (is.null(labels)) {
    labels <- paste0("risk_", seq_len(n_levels))
  }
  if (!is.character(labels) || length(labels) != n_levels || anyNA(labels)) {
    stop(sprintf("`labels` must be NULL or a character vector of %d non-missing names.",
                 n_levels), call. = FALSE)
  }
  if (anyDuplicated(labels)) {
    stop("`labels` must be unique.", call. = FALSE)
  }

  structure(
    list(
      fractions          = as.numeric(fractions),
      relative_risk      = as.numeric(relative_risk),
      reference          = reference,
      labels             = labels,
      n_levels           = n_levels,
      ## Fraction-weighted mean relative risk. This is the factor that turns a
      ## baseline risk into an average per-contact transmission probability, and
      ## hence the factor linking mean contacts to R0.
      mean_relative_risk = sum(fractions * relative_risk),
      ## Largest relative risk, which sets the feasibility ceiling: the top tier's
      ## transmission probability (baseline * max_relative_risk) must stay <= 1.
      max_relative_risk  = max(relative_risk)
    ),
    class = "fiber_contact_risk"
  )
}

#' @export
print.fiber_contact_risk <- function(x, ...) {
  cat("<fiber contact risk structure>\n")
  cat(sprintf("  %d tiers, reference tier: %s\n", x$n_levels, x$labels[x$reference]))
  print(data.frame(
    level         = seq_len(x$n_levels),
    label         = x$labels,
    fraction      = x$fractions,
    relative_risk = x$relative_risk,
    stringsAsFactors = FALSE
  ), row.names = FALSE)
  cat(sprintf("  mean relative risk: %.4f  (R0 = mean contacts x baseline risk x this)\n",
              x$mean_relative_risk))
  cat(sprintf("  max relative risk:  %.4f  (baseline risk must be <= %.4f)\n",
              x$max_relative_risk, 1 / x$max_relative_risk))
  invisible(x)
}

#' Build a graded contact risk structure from a single spread parameter
#'
#' Convenience constructor for a log-spaced gradient of relative risks, so the
#' amount of contact heterogeneity can be swept with one number instead of
#' hand-writing a vector of weights. Relative risks run from 1 (the reference,
#' lowest tier) up to \code{ratio} across \code{n_levels} tiers.
#'
#' @param n_levels Integer, number of tiers. Defaults to 5.
#' @param ratio Positive numeric. Ratio of the top tier's risk to the reference
#'   tier's risk. \code{ratio = 1} reproduces a flat structure.
#' @param fractions Optional numeric vector of tier shares (must sum to 1).
#'   Defaults to equal shares.
#' @param reference Integer index of the reference tier. Defaults to 1 (the
#'   lowest-risk tier).
#' @param labels Optional character vector of tier names.
#'
#' @return An object of class \code{"fiber_contact_risk"}.
#'
#' @examples
#' contact_risk_gradient(n_levels = 5, ratio = 10)
#' @export
contact_risk_gradient <- function(n_levels  = 5,
                                  ratio     = 10,
                                  fractions = NULL,
                                  reference = 1L,
                                  labels    = NULL) {

  if (!is.numeric(n_levels) || length(n_levels) != 1L || is.na(n_levels) ||
      n_levels != round(n_levels) || n_levels < 1) {
    stop("`n_levels` must be a single positive integer.", call. = FALSE)
  }
  n_levels <- as.integer(n_levels)

  if (!is.numeric(ratio) || length(ratio) != 1L || !is.finite(ratio) || ratio <= 0) {
    stop("`ratio` must be a single finite positive numeric value.", call. = FALSE)
  }

  if (is.null(fractions)) {
    fractions <- rep(1 / n_levels, n_levels)
  }

  relative_risk <- if (n_levels == 1L) {
    1
  } else {
    exp(seq(0, log(ratio), length.out = n_levels))
  }

  make_contact_risk(fractions     = fractions,
                    relative_risk = relative_risk,
                    reference     = reference,
                    labels        = labels)
}

## Internal: accept either a fiber_contact_risk object or a plain list of
## make_contact_risk() arguments, and return a validated structure. NULL falls
## back to `default`, which lets each transmission route inherit a shared
## structure unless it overrides it.
as_contact_risk <- function(x, default = NULL, param_name = "contact_risk") {
  if (is.null(x)) {
    if (is.null(default)) {
      return(make_contact_risk())
    }
    return(as_contact_risk(default, NULL, param_name))
  }
  if (inherits(x, "fiber_contact_risk")) {
    return(x)
  }
  if (is.list(x)) {
    return(do.call(make_contact_risk, x))
  }
  stop(sprintf("`%s` must be NULL, a `fiber_contact_risk` object (see make_contact_risk()), or a named list of its arguments.",
               param_name), call. = FALSE)
}

#' Expected secondary infections implied by a contact structure
#'
#' Forward direction of fiber's contact-first calibration:
#' \deqn{R_0 = \bar{c} \times p_0 \times \overline{rr}}
#' where \eqn{\bar{c}} is the mean number of contacts, \eqn{p_0} the baseline
#' per-contact transmission probability (for the reference tier) and
#' \eqn{\overline{rr}} the fraction-weighted mean relative risk.
#'
#' The contact distribution's overdispersion does not appear: thinning a
#' Negative Binomial leaves its mean unchanged, so superspreading can be dialled
#' in the contact distribution without disturbing this calibration.
#'
#' This is the expected offspring count for one transmission route, \emph{before}
#' any of the intervention layers (PPE, hospital quarantine, safe burial) thin it.
#'
#' @param mn_contacts Positive numeric or function(t). Mean of the Negative
#'   Binomial contact distribution.
#' @param baseline_risk Numeric in \code{[0, 1]} or function(t). Per-contact
#'   transmission probability for the reference tier.
#' @param contact_risk A \code{\link{make_contact_risk}} structure (or a named
#'   list of its arguments).
#'
#' @return A single numeric value, or -- if either input is a function of time --
#'   a \code{function(t)} returning the implied R0 at each time.
#'
#' @examples
#' contact_r0_from_baseline_risk(mn_contacts = 20, baseline_risk = 0.1)
#' @export
contact_r0_from_baseline_risk <- function(mn_contacts,
                                          baseline_risk,
                                          contact_risk = make_contact_risk()) {
  risk <- as_contact_risk(contact_risk)

  if (is.function(mn_contacts) || is.function(baseline_risk)) {
    return(function(t) {
      m <- resolve_positive_time_varying(mn_contacts, t, "mn_contacts")
      p <- resolve_time_varying(baseline_risk, t, "baseline_risk")
      m * p * risk$mean_relative_risk
    })
  }

  m <- resolve_positive_time_varying(mn_contacts, 0, "mn_contacts")
  p <- resolve_time_varying(baseline_risk, 0, "baseline_risk")
  m * p * risk$mean_relative_risk
}

#' Solve the baseline per-contact risk from a target R0
#'
#' Inverse of \code{\link{contact_r0_from_baseline_risk}}: given a target
#' expected number of secondary infections, the mean contact number and the risk
#' tier structure, return the baseline per-contact transmission probability that
#' delivers it:
#' \deqn{p_0 = R_0 / (\bar{c} \times \overline{rr})}
#'
#' Every tier's transmission probability must stay in \code{[0, 1]}, so the
#' highest-risk tier caps the achievable R0 at
#' \code{mn_contacts * mean_relative_risk / max_relative_risk} (see
#' \code{\link{max_contact_r0}}). Asking for more than that is an error rather
#' than a silently clipped top tier.
#'
#' @param r0 Non-negative numeric or function(t). Target expected secondary
#'   infections for this transmission route, before intervention thinning.
#' @param mn_contacts Positive numeric or function(t). Mean of the Negative
#'   Binomial contact distribution.
#' @param contact_risk A \code{\link{make_contact_risk}} structure (or a named
#'   list of its arguments).
#' @param param_name Character, used to make error messages point at the
#'   offending simulation argument.
#'
#' @return A single numeric baseline risk in \code{[0, 1]}, or -- if either input
#'   is a function of time -- a \code{function(t)} returning it at each time.
#'
#' @examples
#' ## 20 contacts on average, want R0 = 2 with a 10-fold risk gradient
#' risk <- contact_risk_gradient(n_levels = 5, ratio = 10)
#' baseline_risk_from_contact_r0(r0 = 2, mn_contacts = 20, contact_risk = risk)
#' @export
baseline_risk_from_contact_r0 <- function(r0,
                                          mn_contacts,
                                          contact_risk = make_contact_risk(),
                                          param_name   = "r0") {
  risk <- as_contact_risk(contact_risk)

  solve_one <- function(r0_t, m_t, t_label = NULL) {
    if (any(!is.finite(r0_t)) || any(r0_t < 0)) {
      stop(sprintf("`%s` must resolve to finite, non-negative value(s).", param_name),
           call. = FALSE)
    }
    p0 <- r0_t / (m_t * risk$mean_relative_risk)
    ceiling_r0 <- m_t * risk$mean_relative_risk / risk$max_relative_risk
    bad <- p0 * risk$max_relative_risk > 1 + 1e-12
    if (any(bad)) {
      i <- which(bad)[1]
      stop(sprintf(
        paste0("`%s` = %s is not achievable with mean contacts = %s and this risk structure: ",
               "it needs a baseline risk of %s, which puts the highest-risk tier at %s (must be <= 1). ",
               "The largest achievable R0 here is %s -- raise the mean contact number, ",
               "or narrow the relative-risk spread.%s"),
        param_name,
        format(round(r0_t[i], 4)),
        format(round(m_t[i], 4)),
        format(round(p0[i], 6)),
        format(round(p0[i] * risk$max_relative_risk, 4)),
        format(round(ceiling_r0[i], 4)),
        if (is.null(t_label)) "" else sprintf(" (first failure at t = %s)", format(round(t_label[i], 3)))
      ), call. = FALSE)
    }
    p0
  }

  if (is.function(r0) || is.function(mn_contacts)) {
    return(function(t) {
      r0_t <- resolve_time_varying(r0, t, param_name)
      m_t  <- resolve_positive_time_varying(mn_contacts, t, "mn_contacts")
      solve_one(r0_t, m_t, t)
    })
  }

  m <- resolve_positive_time_varying(mn_contacts, 0, "mn_contacts")
  r <- resolve_time_varying(r0, 0, param_name)
  solve_one(r, m)
}

#' Largest R0 a contact structure can produce
#'
#' The feasibility ceiling implied by \code{\link{baseline_risk_from_contact_r0}}:
#' the highest-risk tier's transmission probability cannot exceed 1, so
#' \deqn{R_0^{max} = \bar{c} \times \overline{rr} / \max(rr)}
#' With a flat risk structure this is just the mean contact number -- you cannot
#' infect more people than you meet.
#'
#' @param mn_contacts Positive numeric. Mean number of contacts.
#' @param contact_risk A \code{\link{make_contact_risk}} structure (or a named
#'   list of its arguments).
#'
#' @return A single numeric value.
#' @export
max_contact_r0 <- function(mn_contacts, contact_risk = make_contact_risk()) {
  risk <- as_contact_risk(contact_risk)
  m <- resolve_positive_time_varying(mn_contacts, 0, "mn_contacts")
  m * risk$mean_relative_risk / risk$max_relative_risk
}

#' Rough overall R0 for an index case, including the funeral route
#'
#' A back-of-envelope assembly of the per-route expectations into a single
#' number for an index general-population case: the direct (community and
#' hospital) route, plus the funeral route weighted by the chance the case dies
#' and receives an unsafe burial.
#'
#' \strong{This is a parameter-set sanity check, not an estimate of the simulated
#' R.} It is deliberately "rough":
#' \itemize{
#'   \item It is \emph{unmitigated} -- computed before PPE, hospital quarantine
#'         and safe-burial thinning, all of which reduce the realised R.
#'   \item It treats the probability of dying crudely as
#'         \code{prob_symptomatic * prob_death_comm}, whereas the simulation
#'         resolves death jointly with whether the case reaches hospital in time.
#'   \item It ignores the healthcare-worker parent route, which applies only to
#'         cases who are themselves HCWs.
#' }
#' Use \code{\link{compute_reproduction_number}} on a simulated tree for the
#' realised reproduction number.
#'
#' @param mn_contacts_genPop,baseline_risk_genPop Contact mean and baseline risk
#'   for the general-population transmission route.
#' @param mn_contacts_funeral,baseline_risk_funeral Contact mean and baseline
#'   risk for the funeral route. Both \code{NULL} drops the funeral term.
#' @param contact_risk Risk structure shared by both routes unless overridden.
#' @param contact_risk_genPop,contact_risk_funeral Optional per-route risk
#'   structures; \code{NULL} inherits \code{contact_risk}.
#' @param prob_symptomatic,prob_death_comm Probabilities used for the crude
#'   death weighting of the funeral term.
#' @param p_unsafe_funeral Numeric in \code{[0, 1]} or function(t). Probability a
#'   death leads to an unsafe funeral.
#' @param t Numeric time at which to evaluate any time-varying inputs. Defaults
#'   to 0.
#'
#' @return A list with the \code{direct} and \code{funeral} contributions, their
#'   \code{total}, and the inputs used (\code{p_death}, \code{p_unsafe_funeral},
#'   \code{t}).
#' @export
compute_rough_r0 <- function(mn_contacts_genPop,
                             baseline_risk_genPop,
                             mn_contacts_funeral  = NULL,
                             baseline_risk_funeral = NULL,
                             contact_risk          = make_contact_risk(),
                             contact_risk_genPop   = NULL,
                             contact_risk_funeral  = NULL,
                             prob_symptomatic      = 1,
                             prob_death_comm       = 0,
                             p_unsafe_funeral      = 0,
                             t                     = 0) {

  base_risk    <- as_contact_risk(contact_risk)
  risk_genPop  <- as_contact_risk(contact_risk_genPop,  base_risk, "contact_risk_genPop")
  risk_funeral <- as_contact_risk(contact_risk_funeral, base_risk, "contact_risk_funeral")

  at <- function(param, name, positive = FALSE) {
    if (positive) resolve_positive_time_varying(param, t, name)
    else          resolve_time_varying(param, t, name)
  }

  direct <- at(mn_contacts_genPop, "mn_contacts_genPop", positive = TRUE) *
            at(baseline_risk_genPop, "baseline_risk_genPop") *
            risk_genPop$mean_relative_risk

  p_death <- at(prob_symptomatic, "prob_symptomatic") * at(prob_death_comm, "prob_death_comm")
  p_unsafe <- at(p_unsafe_funeral, "p_unsafe_funeral")

  funeral <- if (is.null(mn_contacts_funeral) || is.null(baseline_risk_funeral)) {
    0
  } else {
    p_death * p_unsafe *
      at(mn_contacts_funeral, "mn_contacts_funeral", positive = TRUE) *
      at(baseline_risk_funeral, "baseline_risk_funeral") *
      risk_funeral$mean_relative_risk
  }

  list(
    direct           = direct,
    funeral          = funeral,
    total            = direct + funeral,
    p_death          = p_death,
    p_unsafe_funeral = p_unsafe,
    t                = t
  )
}

## ---------------------------------------------------------------------------
## Runtime pieces used by the offspring functions
## ---------------------------------------------------------------------------

## Decide whether an offspring function runs in contact mode or legacy offspring
## mode, and validate the parameters that mode needs. Exactly one of the two
## means must be supplied. `suffix` names the route ("genPop", "hcw",
## "funeral") so error messages point at the actual argument.
resolve_transmission_mode <- function(mn_offspring,
                                      overdisp_offspring,
                                      mn_contacts,
                                      overdisp_contacts,
                                      baseline_risk,
                                      contact_risk,
                                      suffix) {

  nm <- function(stem) paste0(stem, "_", suffix)

  positive_or_time_varying <- function(param, param_name) {
    if (is.function(param)) return(invisible(NULL))
    if (!is.numeric(param) || length(param) != 1L || is.na(param) || param <= 0) {
      stop(sprintf("`%s` must be a function(t) or a single positive numeric value.", param_name),
           call. = FALSE)
    }
    invisible(NULL)
  }

  positive_scalar <- function(param, param_name) {
    if (is.null(param) || !is.numeric(param) || length(param) != 1L ||
        is.na(param) || param <= 0) {
      stop(sprintf("`%s` must be a single positive numeric value.", param_name), call. = FALSE)
    }
    invisible(NULL)
  }

  has_offspring <- !is.null(mn_offspring)
  has_contacts  <- !is.null(mn_contacts)

  if (has_offspring && has_contacts) {
    stop(sprintf("Supply either `%s` (offspring mode) or `%s` (contact mode), not both.",
                 nm("mn_offspring"), nm("mn_contacts")), call. = FALSE)
  }
  if (!has_offspring && !has_contacts) {
    stop(sprintf("Supply one of `%s` (offspring mode) or `%s` (contact mode).",
                 nm("mn_offspring"), nm("mn_contacts")), call. = FALSE)
  }

  if (has_contacts) {
    positive_or_time_varying(mn_contacts, nm("mn_contacts"))
    positive_scalar(overdisp_contacts, nm("overdisp_contacts"))
    if (is.null(baseline_risk)) {
      stop(sprintf("`%s` is required in contact mode (solve it from a target R0 with baseline_risk_from_contact_r0()).",
                   nm("baseline_risk")), call. = FALSE)
    }
    if (!is.function(baseline_risk) &&
        (!is.numeric(baseline_risk) || length(baseline_risk) != 1L ||
         is.na(baseline_risk) || baseline_risk < 0 || baseline_risk > 1)) {
      stop(sprintf("`%s` must be a function(t) or a single numeric in [0, 1].", nm("baseline_risk")),
           call. = FALSE)
    }
    return(list(
      mode          = "contact",
      mn            = mn_contacts,
      size          = overdisp_contacts,
      baseline_risk = baseline_risk,
      risk          = as_contact_risk(contact_risk, NULL, nm("contact_risk"))
    ))
  }

  positive_or_time_varying(mn_offspring, nm("mn_offspring"))
  positive_scalar(overdisp_offspring, nm("overdisp_offspring"))
  list(mode = "offspring", mn = mn_offspring, size = overdisp_offspring)
}

## Assign a risk tier to each of `n` contacts, drawn independently from the
## structure's tier fractions. Consumes exactly one RNG call.
draw_contact_risk_levels <- function(n, contact_risk) {
  if (n == 0L) return(integer(0))
  sample.int(contact_risk$n_levels, size = n, replace = TRUE,
             prob = contact_risk$fractions)
}

## Per-contact transmission probability = baseline_risk(t) * relative_risk[tier].
## `baseline_risk` may be a scalar or a function of absolute calendar time, and
## is resolved at each contact's own contact time (matching how PPE coverage and
## quarantine efficacy are resolved per event rather than per parent).
contact_transmission_prob <- function(levels,
                                      contact_times_absolute,
                                      baseline_risk,
                                      contact_risk,
                                      param_name = "baseline_risk") {
  if (length(levels) == 0L) return(numeric(0))
  p0 <- resolve_time_varying(baseline_risk, contact_times_absolute, param_name)
  if (any(p0 < 0 | p0 > 1)) {
    stop(sprintf("`%s` must resolve to value(s) in [0, 1].", param_name), call. = FALSE)
  }
  p <- p0 * contact_risk$relative_risk[levels]
  if (any(p > 1 + 1e-12)) {
    i <- which(p > 1 + 1e-12)[1]
    stop(sprintf(
      "`%s` = %s at t = %s gives tier '%s' (relative risk %s) a transmission probability of %s, which exceeds 1. Lower the baseline risk or narrow the relative-risk spread.",
      param_name,
      format(round(p0[i], 6)),
      format(round(contact_times_absolute[i], 3)),
      contact_risk$labels[levels[i]],
      format(round(contact_risk$relative_risk[levels[i]], 4)),
      format(round(p[i], 4))
    ), call. = FALSE)
  }
  pmin(p, 1)
}

## Empty (0-row) contact log, matching the columns built by new_contact_log().
empty_contact_log <- function() {
  data.frame(
    parent                  = integer(0),
    case_id                 = integer(0),
    record_type             = character(0),
    class                   = character(0),
    infection_location      = character(0),
    time_contact_relative   = numeric(0),
    time_contact_absolute   = numeric(0),
    risk_level              = integer(0),
    risk_category           = character(0),
    relative_risk           = numeric(0),
    transmission_prob       = numeric(0),
    blocked_by              = character(0),
    stringsAsFactors        = FALSE
  )
}

## Build the per-parent contact log. One row per contact generated, in the order
## the contacts were drawn.
##
##   record_type  "infection" if the contact ended up as a realised case,
##                "contact" otherwise.
##   blocked_by   NA for realised infections; otherwise the layer that stopped
##                it -- "no_transmission" (the risk draw), the route's
##                intervention layer, or "obv_pep".
##
## `case_id` is left NA here because ids are assigned by the simulation loop;
## branching_process_main fills it in for the "infection" rows, which appear in
## the same order as the offspring rows this call returns.
new_contact_log <- function(parent_id,
                            offspring_class,
                            infection_location,
                            time_contact_relative,
                            time_contact_absolute,
                            risk_levels,
                            transmission_prob,
                            contact_risk,
                            transmitted,
                            kept,
                            realised,
                            intervention_label) {

  n <- length(offspring_class)
  if (n == 0L) return(empty_contact_log())

  blocked_by <- rep(NA_character_, n)
  blocked_by[!transmitted]           <- "no_transmission"
  blocked_by[transmitted & !kept]    <- intervention_label
  blocked_by[kept & !realised]       <- "obv_pep"

  data.frame(
    parent                = rep(parent_id, n),
    case_id               = rep(NA_integer_, n),
    record_type           = ifelse(realised, "infection", "contact"),
    class                 = offspring_class,
    infection_location    = infection_location,
    time_contact_relative = time_contact_relative,
    time_contact_absolute = time_contact_absolute,
    risk_level            = as.integer(risk_levels),
    risk_category         = contact_risk$labels[risk_levels],
    relative_risk         = contact_risk$relative_risk[risk_levels],
    transmission_prob     = transmission_prob,
    blocked_by            = blocked_by,
    stringsAsFactors      = FALSE
  )
}
