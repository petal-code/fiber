## Ring vaccination for the fiber branching-process model.
##
## Ring vaccination is a TARGETED, CROSS-GENERATIONAL intervention, implemented as
## two main-loop hooks around the existing offspring-generation step (see
## branching_process_main):
##
##   Hook A (parent infectiousness): a vaccinated-breakthrough parent transmits
##     less, so its realised brood is thinned by (1 - efficacy_transmission).
##   Hook B (protect the children): apply_ring_vax_gate() decides, for each
##     candidate child the parent generated, whether ring vaccination prevents
##     that infection (averting it before complete_offspring_info() resolves any
##     natural history), and records vaccination / breakthrough status.
##
## The model is parent-centric rather than index-centric (see the design notes in
## the package docs): every parent, when expanded, vaccinates its own children
## using the campaigns it carries -- its OWN campaign (ring 1, if it is a detected
## index) plus its PARENT's campaign (ring 2, the grandparent of these children).
## Because fiber expands infections in strict infection-time order, both campaigns
## are fully known the moment a parent is expanded, so each child's fate is decided
## once, at birth, with no forward simulation.
##
## CLOCKS. Three distinct times drive a campaign anchored at case A:
##   * probability clock  -- A's INFECTION time: trace_prob, coverage, detection
##                           and the delays are resolved here (campaign intensity
##                           is attributed to the triggering case's infection time).
##   * logistics clock    -- A's DETECTION time (= symptom onset + reporting_delay):
##                           vaccination becomes available at v = detection + logistical_delay
##                           (+ a per-hop increment for ring 2), protection at p = v + protection_delay.
##   * race clock         -- each CHILD's infection time, compared against v (reached
##                           in time?) and p (protected in time?).
##
## CHAIN, per ring: traced ~ Bernoulli(trace_prob) -> received ~ Bernoulli(coverage | traced)
##   -> protected-in-time iff p <= child infection time -> averted ~ Bernoulli(efficacy_infection).
## A child is AVERTED if any ring averts it; a retained child that was protected in
## time (efficacy failed) is a BREAKTHROUGH (reduced transmitter, Hook A next
## generation); a retained child vaccinated but infected before protection is a full
## transmitter.
##
## "Require intermediate traced" (ring 2): the grandparent's campaign reaches a
## grandchild only if the intermediate parent was itself traced by the grandparent's
## ring 1. Tracing is separate from vaccination, so an intermediate found-but-unvaccinated
## still relays the ring.

## Per-call ring-vaccination counters, summed across offspring-generation steps in
## branching_process_main and surfaced by summarise_output(). All four gate counters
## are UNIQUE-CONTACT counts among the parent's realised would-be offspring (a
## lower-bound "dose proxy": real campaigns also vaccinate contacts who would never
## have been infected, which a branching process without a susceptible pool cannot
## represent). They are nested per call: prevented <= protected <= vaccinated <= traced.
##
## `prevented_deaths` is NOT a gate counter: like the OBV analogue it is the subset
## of `prevented` infections that WOULD have died had they occurred, resolved once
## after the simulation loop by replaying the averted infections through the same
## outcome model as realised cases (see branching_process_main). It stays 0L here.
empty_ring_vax_num_treated <- function() {
  list(
    traced           = 0L,   # contacts traced by >= 1 ring (and would have been infected)
    vaccinated       = 0L,   # contacts who received >= 1 dose (traced AND covered)
    protected        = 0L,   # vaccinated contacts protected in time (breakthrough-eligible)
    prevented        = 0L,   # contacts whose infection was averted (protected AND efficacy took)
    prevented_deaths = 0L    # subset of `prevented` that would have died (deferred counterfactual)
  )
}

## Empty (0-row) container for the per-call snapshot of infections ring vaccination
## prevented. Carries just enough to resolve each one's counterfactual death later:
## class (drives hospitalisation/CFR), infection location, and the absolute infection
## time (the calendar clock for time-varying parameters). Mirrors the OBV analogue.
empty_ring_vax_prevented_info <- function() {
  data.frame(
    class                   = character(0),
    infection_location      = character(0),
    time_infection_absolute = numeric(0),
    stringsAsFactors        = FALSE
  )
}

#' Apply the ring-vaccination gate to a parent's candidate brood (Hook B)
#'
#' Decides, for every candidate child a parent generated, whether ring vaccination
#' prevents that infection, using the parent's own campaign (ring 1) and the
#' grandparent's campaign (ring 2). All RNG drawn here is consumed only when ring
#' vaccination is enabled (the caller guards the call), so a disabled run is
#' bit-for-bit identical to the pre-feature simulation.
#'
#' @param brood Data frame of candidate offspring (post PPE/ETU/OBV thinning and
#'   post Hook-A breakthrough thinning). Reads `class`, `infection_location` and
#'   `time_infection_relative` (relative to the parent's infection time).
#' @param parent_time_infection_absolute Scalar. The parent's absolute infection
#'   time (ring-1 probability clock; also converts relative child times to absolute).
#' @param parent_symptomatic Logical scalar. Whether the parent is symptomatic
#'   (only symptomatic parents can be detected and start a ring-1 campaign).
#' @param parent_time_symptom_onset_absolute Scalar (NA if asymptomatic). The
#'   parent's absolute symptom-onset time (the ring-1 logistics clock origin).
#' @param grandparent_detection_time Scalar or NA. The grandparent's detection
#'   time (ring-2 logistics clock origin); NA when there is no detected grandparent.
#' @param grandparent_infection_time Scalar or NA. The grandparent's infection time
#'   (ring-2 probability clock).
#' @param parent_traced_by_grandparent Logical scalar. Whether the parent was traced
#'   by the grandparent's ring-1 campaign; gates ring 2 when
#'   `ring_vax_require_intermediate_traced = TRUE`.
#' @param ring_vax_start Scalar. Calendar time from which vaccination is available;
#'   a ring acts on a child only if the vaccination time is at or after this.
#' @param ring_vax_n_rings Integer 1 or 2. Number of rings (1 disables ring 2).
#' @param ring_vax_detection_prob,ring_vax_reporting_delay,ring_vax_trace_prob,ring_vax_coverage,ring_vax_logistical_delay,ring_vax_protection_delay,ring_vax_ring2_delay_increment
#'   Scalars or `function(t)`; resolved on the probability clock described above.
#' @param ring_vax_efficacy_infection Scalar in [0, 1]. All-or-nothing vaccine
#'   efficacy against infection (applied once per protected-in-time contact).
#' @param ring_vax_independent_coverage Logical. TRUE: an independent coverage draw
#'   per ring (a contact missed by one ring can still be vaccinated by the other).
#'   FALSE: one acceptance draw per contact, shared across rings.
#' @param ring_vax_require_intermediate_traced Logical. If TRUE, ring 2 reaches a
#'   grandchild only when the intermediate parent was traced by the grandparent.
#' @param ring_vax_target_class,ring_vax_target_locations Character vectors of
#'   eligible offspring classes / infection locations.
#'
#' @return A list with per-child vectors (length `nrow(brood)`): `keep`
#'   (FALSE = averted), `vaccinated`, `breakthrough`, `time_vaccinated`,
#'   `traced_by_parent` (ring-1 trace outcome, stamped on surviving children); plus
#'   scalars `parent_detection_time` (this parent's detection time, stamped on
#'   surviving children as their ring-2 grandparent anchor), the `num_treated`
#'   counter list, and `prevented_info` (a frame of averted infections for the
#'   deferred prevented-deaths counterfactual).
#' @noRd
apply_ring_vax_gate <- function(brood,
                                parent_time_infection_absolute,
                                parent_symptomatic,
                                parent_time_symptom_onset_absolute,
                                grandparent_detection_time = NA_real_,
                                grandparent_infection_time = NA_real_,
                                parent_traced_by_grandparent = FALSE,
                                ring_vax_start = 0,
                                ring_vax_n_rings = 2,
                                ring_vax_detection_prob = 1,
                                ring_vax_reporting_delay = 0,
                                ring_vax_trace_prob = 0,
                                ring_vax_coverage = 0,
                                ring_vax_efficacy_infection = 0,
                                ring_vax_logistical_delay = 0,
                                ring_vax_protection_delay = 0,
                                ring_vax_ring2_delay_increment = 0,
                                ring_vax_independent_coverage = TRUE,
                                ring_vax_require_intermediate_traced = TRUE,
                                ring_vax_target_class = c("genPop", "HCW"),
                                ring_vax_target_locations = c("community", "hospital", "funeral")) {

  n <- nrow(brood)

  ## Per-child outputs (default = untouched / not vaccinated / retained).
  keep             <- rep(TRUE,  n)
  vaccinated       <- rep(FALSE, n)
  breakthrough     <- rep(FALSE, n)
  time_vaccinated  <- rep(NA_real_, n)
  traced_by_parent <- rep(FALSE, n)   # ring-1 trace outcome -> stamped on surviving children

  num_treated    <- empty_ring_vax_num_treated()
  prevented_info <- empty_ring_vax_prevented_info()

  ## --- Parent detection (the ring-1 anchor) --------------------------------
  ## The parent becomes a detected index iff it is symptomatic AND a detection
  ## draw succeeds. Detection probability and reporting delay are resolved on the
  ## parent's infection clock. Detection time = symptom onset + reporting delay.
  ## This is the FIRST RNG draw in the gate (a single Bernoulli) so the draw order
  ## is well-defined for reproducibility.
  parent_detection_time <- NA_real_
  if (isTRUE(parent_symptomatic) && !is.na(parent_time_symptom_onset_absolute)) {
    det_prob <- resolve_time_varying(ring_vax_detection_prob,
                                     parent_time_infection_absolute,
                                     "ring_vax_detection_prob")
    if (rbinom(n = 1, size = 1, prob = det_prob) == 1L) {
      rep_delay <- resolve_time_varying(ring_vax_reporting_delay,
                                        parent_time_infection_absolute,
                                        "ring_vax_reporting_delay")
      parent_detection_time <- parent_time_symptom_onset_absolute + rep_delay
    }
  }

  empty_result <- function() {
    list(keep = keep, vaccinated = vaccinated, breakthrough = breakthrough,
         time_vaccinated = time_vaccinated, traced_by_parent = traced_by_parent,
         parent_detection_time = parent_detection_time,
         num_treated = num_treated, prevented_info = prevented_info)
  }
  if (n == 0L) return(empty_result())

  t_inf_child <- parent_time_infection_absolute + brood$time_infection_relative
  eligible    <- (brood$class %in% ring_vax_target_class) &
                 (brood$infection_location %in% ring_vax_target_locations)

  ## --- Ring 1: anchor = the parent -----------------------------------------
  ## Active iff the parent is a detected index. Vaccination time v1 = detection +
  ## logistical delay; protection time p1 = v1 + protection delay (delays on the
  ## parent's infection clock). A child is "reached" if it is eligible, the
  ## campaign has started (v1 >= ring_vax_start) and it is vaccinated before it is
  ## infected (v1 <= child infection time).
  traced1 <- rep(FALSE, n)
  v1 <- NA_real_; p1 <- NA_real_
  if (!is.na(parent_detection_time)) {
    ld1 <- resolve_time_varying(ring_vax_logistical_delay,
                                parent_time_infection_absolute, "ring_vax_logistical_delay")
    pd1 <- resolve_time_varying(ring_vax_protection_delay,
                                parent_time_infection_absolute, "ring_vax_protection_delay")
    v1 <- parent_detection_time + ld1
    p1 <- v1 + pd1
    if (v1 >= ring_vax_start) {
      reach1 <- eligible & (v1 <= t_inf_child)
      idx1 <- which(reach1)
      if (length(idx1) > 0L) {
        tp1 <- resolve_time_varying(ring_vax_trace_prob,
                                    parent_time_infection_absolute, "ring_vax_trace_prob")
        traced1[idx1] <- rbinom(n = length(idx1), size = 1, prob = tp1) == 1L
      }
    }
  }
  ## Surviving children inherit their ring-1 trace outcome: it gates whether THIS
  ## parent's campaign can later reach their children as ring 2.
  traced_by_parent <- traced1

  ## --- Ring 2: anchor = the grandparent ------------------------------------
  ## Active iff a second ring is requested, a detected grandparent exists, and
  ## (when required) the parent was itself traced by the grandparent. Timing uses
  ## the grandparent's detection time plus an optional second-hop increment;
  ## probabilities/delays are resolved on the grandparent's infection clock.
  traced2 <- rep(FALSE, n)
  v2 <- NA_real_; p2 <- NA_real_
  ring2_active <- (ring_vax_n_rings >= 2) &&
                  !is.na(grandparent_detection_time) &&
                  (!isTRUE(ring_vax_require_intermediate_traced) ||
                     isTRUE(parent_traced_by_grandparent))
  if (ring2_active) {
    ld2 <- resolve_time_varying(ring_vax_logistical_delay,
                                grandparent_infection_time, "ring_vax_logistical_delay")
    pd2 <- resolve_time_varying(ring_vax_protection_delay,
                                grandparent_infection_time, "ring_vax_protection_delay")
    inc <- resolve_time_varying(ring_vax_ring2_delay_increment,
                                grandparent_infection_time, "ring_vax_ring2_delay_increment")
    v2 <- grandparent_detection_time + ld2 + inc
    p2 <- v2 + pd2
    if (v2 >= ring_vax_start) {
      reach2 <- eligible & (v2 <= t_inf_child)
      idx2 <- which(reach2)
      if (length(idx2) > 0L) {
        tp2 <- resolve_time_varying(ring_vax_trace_prob,
                                    grandparent_infection_time, "ring_vax_trace_prob")
        traced2[idx2] <- rbinom(n = length(idx2), size = 1, prob = tp2) == 1L
      }
    }
  }

  ## --- Coverage (vaccine receipt | traced) ---------------------------------
  ## Independent: an independent coverage draw per ring (two chances to receive).
  ## Shared: one acceptance draw per contact (resolved on the parent's clock),
  ## reused across rings -- a refuser refuses every ring.
  vacc1 <- rep(FALSE, n); vacc2 <- rep(FALSE, n)
  if (isTRUE(ring_vax_independent_coverage)) {
    tr1 <- which(traced1)
    if (length(tr1) > 0L) {
      cov1 <- resolve_time_varying(ring_vax_coverage,
                                   parent_time_infection_absolute, "ring_vax_coverage")
      vacc1[tr1] <- rbinom(n = length(tr1), size = 1, prob = cov1) == 1L
    }
    tr2 <- which(traced2)
    if (length(tr2) > 0L) {
      cov2 <- resolve_time_varying(ring_vax_coverage,
                                   grandparent_infection_time, "ring_vax_coverage")
      vacc2[tr2] <- rbinom(n = length(tr2), size = 1, prob = cov2) == 1L
    }
  } else {
    traced_any <- traced1 | traced2
    ta <- which(traced_any)
    accept <- rep(FALSE, n)
    if (length(ta) > 0L) {
      covp <- resolve_time_varying(ring_vax_coverage,
                                   parent_time_infection_absolute, "ring_vax_coverage")
      accept[ta] <- rbinom(n = length(ta), size = 1, prob = covp) == 1L
    }
    vacc1 <- traced1 & accept
    vacc2 <- traced2 & accept
  }

  vaccinated <- vacc1 | vacc2

  ## Effective dose / protection times = earliest across the rings that vaccinated
  ## the contact (v1, v2, p1, p2 are scalars; an inactive ring never vaccinates).
  v_eff <- rep(Inf, n); p_eff <- rep(Inf, n)
  if (any(vacc1)) { v_eff[vacc1] <- pmin(v_eff[vacc1], v1); p_eff[vacc1] <- pmin(p_eff[vacc1], p1) }
  if (any(vacc2)) { v_eff[vacc2] <- pmin(v_eff[vacc2], v2); p_eff[vacc2] <- pmin(p_eff[vacc2], p2) }
  time_vaccinated[vaccinated] <- v_eff[vaccinated]

  protected_in_time <- vaccinated & (p_eff <= t_inf_child)

  ## --- Efficacy (one all-or-nothing draw per protected-in-time contact) -----
  prot_idx <- which(protected_in_time)
  if (length(prot_idx) > 0L) {
    averted <- rbinom(n = length(prot_idx), size = 1, prob = ring_vax_efficacy_infection) == 1L
    keep[prot_idx[averted]]          <- FALSE
    breakthrough[prot_idx[!averted]] <- TRUE
  }

  ## --- Counters and the averted-infection snapshot --------------------------
  num_treated$traced     <- sum(traced1 | traced2)
  num_treated$vaccinated <- sum(vaccinated)
  num_treated$protected  <- sum(protected_in_time)
  num_treated$prevented  <- sum(!keep)

  av_idx <- which(!keep)
  if (length(av_idx) > 0L) {
    prevented_info <- data.frame(
      class                   = brood$class[av_idx],
      infection_location      = brood$infection_location[av_idx],
      time_infection_absolute = t_inf_child[av_idx],
      stringsAsFactors        = FALSE
    )
  }

  list(keep = keep, vaccinated = vaccinated, breakthrough = breakthrough,
       time_vaccinated = time_vaccinated, traced_by_parent = traced_by_parent,
       parent_detection_time = parent_detection_time,
       num_treated = num_treated, prevented_info = prevented_info)
}
