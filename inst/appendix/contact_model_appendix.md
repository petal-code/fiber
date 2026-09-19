# Supplementary Appendix

## Contact-structured transmission, risk-tiered exposure, and contact tracing in a branching-process model of filovirus outbreak control

### Table of Contents

**1. Supplementary Methods**

- Overview of model revisions
- Contact generation and risk-structured transmission
- Contact risk tiers and their parameterisation
- Reduction to the unstructured model, and what heterogeneity does and does not add
- The risk-tier composition of cases differs from that of contacts
- Contact tracing
- Effect of tracing on the probability and timing of hospitalisation
- Calibration of the basic reproduction number
- Presymptomatic transmission
- Model outputs and run diagnostics

**2. Supplementary Tables**

- Supplementary Table S1. Parameters introduced or redefined by the contact structure
- Supplementary Table S2. Illustrative risk-tier structure and its derived quantities

**3. Known limitations of the current implementation**

---

## 1. Supplementary Methods

### Overview of model revisions

The model previously generated secondary infections directly: for each infectious individual, a Negative Binomial draw returned the number of people that individual infected, and each intervention acted by deleting infections from that draw. Exposure events that did not result in transmission were never represented. This was adequate for interventions that act on an infectious case, such as hospital quarantine or safe burial, but it provided no denominator for interventions that act on *exposed contacts*, such as contact tracing or post-exposure prophylaxis, and it folded all heterogeneity in exposure intensity invisibly into a single offspring mean.

We restructured transmission so that the Negative Binomial draw returns the number of **contacts** — exposure events — and infections are then sampled from those contacts. Each contact is assigned to one of a small number of risk tiers, and the tier determines both the probability that the contact results in transmission and the probability that the contact is successfully traced. Contacts that transmit pass unchanged through the model's existing intervention layers (personal protective equipment and infection prevention and control for health workers, post-admission hospital or Ebola treatment unit quarantine, safe burial practices, and post-exposure prophylaxis).

The revision is deliberately constructed so that it reduces exactly to the previous model under a flat risk structure, which we establish formally below. All quantities are defined per transmission route — general population, health-care worker, and funeral — and are parameterised independently for each.

### Contact generation and risk-structured transmission

For an infectious individual *i* of route *R*, the number of contacts is drawn as

```
    N_i ~ NegBin(mean = m_R(t_i), size = k_R)
```

where `m_R` is the mean number of contacts generated over the individual's infectious period and `k_R` is the overdispersion parameter. The mean may be specified as a scalar or as a function of absolute calendar time, resolved once per individual immediately before the draw: at the individual's infection time for the general-population and health-care-worker routes, and at their death time for the funeral route, since the funeral occurs then.

The timing of each contact is drawn from the route's Gamma generation-time distribution, truncated to the interval between the start of infectiousness and the individual's outcome (death or recovery):

```
    t_ij ~ TruncGamma(shape_R, rate_R; lower = a_i, upper = T_i^outcome),   j = 1 … N_i
```

with `a_i = 0` by default and `a_i` equal to the individual's incubation period when presymptomatic transmission is disabled (see below). Truncation is performed exactly, by inverse-CDF sampling on the interval, which is distributionally identical to rejection sampling with unlimited retries.

Each contact *j* is independently allocated to a risk tier *l* ∈ {1 … L} with probability *f_l*, and transmits with probability

```
    P(transmission | tier l, time t) = p_R(t) · r_l
```

where `p_R(t)` is the route's baseline per-contact transmission probability, itself permitted to vary with calendar time, and `r_l` is the tier's relative risk. Critically, `p_R` is resolved at **that contact's own** calendar time rather than at the infector's infection time, so a time-varying transmissibility acts on exposures as they occur.

The setting of each contact (community, hospital, or funeral) and the class of the person exposed (general population or health-care worker) are determined as in the previous model, and the resulting infections then pass through the intervention layers unchanged.

### Contact risk tiers and their parameterisation

A risk structure comprises the tier fractions `f = (f_1 … f_L)` with Σ f_l = 1, the relative risks `r = (r_1 … r_L)`, a per-tier tracing probability `τ = (τ_1 … τ_L)`, and a nominated reference tier. The default is five equally sized tiers with flat relative risk and no tracing; any number of tiers is permitted. A convenience constructor generates a log-spaced gradient from a single spread parameter, optionally with tracing probability interpolated across tiers.

Relative risks are normalised so that the reference tier has relative risk exactly 1. **The reference is the highest-risk tier by default**, so that relative risks fall in (0, 1] and the baseline risk `p_R` is interpretable as the per-contact transmission probability of the riskiest exposure — a household or caregiving contact, which is the quantity that reported secondary attack rates measure. Three considerations motivated this convention:

First, the feasibility requirement that no tier's transmission probability exceed one, `p_R · max_l r_l ≤ 1`, reduces to `p_R ≤ 1` and can therefore never bind. Anchoring instead on the lowest-risk tier leaves an upper bound on `p_R` that moves with the relative-risk spread; with risks spanning a 25-fold range, for example, the baseline would be silently capped at 0.04.

Second, `p_R` and the relative risks then occupy [0, 1] independently, which gives approximate Bayesian computation a rectangular parameter space rather than one whose bounds shift as the spread is itself fitted.

Third, the baseline risk is anchored to a measurable quantity rather than to the rate of transmission at casual contact, which is not directly observable.

Nominating a different reference tier is permitted and constitutes a pure reparameterisation: only the product `p_R · r_l` ever enters the simulation, so rescaling both cannot alter a simulated outcome. This invariance is verified by regression test. The feasibility constraint is evaluated across the simulation horizon before the run begins rather than enforced by clipping during it, so an infeasible specification fails immediately and visibly.

### Reduction to the unstructured model, and what heterogeneity does and does not add

Because tiers are allocated independently across contacts, the *marginal* transmission probability of a contact, before its tier is known, is

```
    p_R · r̄,    where    r̄ = Σ_l f_l r_l
```

irrespective of which tier it lands in. Independent Bernoulli thinning of a Negative Binomial variate returns a Negative Binomial variate with the same overdispersion, so the number of infections generated by an individual is

```
    N_i^inf ~ NegBin(mean = m_R · p_R · r̄,  size = k_R)
```

exactly. Two consequences follow, both of which matter for interpretation.

Under a flat tier structure with `p_R = 1`, the model reduces *identically* to the previous direct offspring draw with `m_R` equal to the former offspring mean. This is not an approximation, and it is the reason every pre-existing test suite passes after the revision with only a parameter rename applied and no assertion altered.

Independent per-contact heterogeneity contributes no additional variance in offspring number and therefore generates no superspreading of its own; overdispersion remains governed entirely by `k_R`. The risk tiers acquire epidemiological consequence only once something downstream conditions on them — which, in this model, is contact tracing.

### The risk-tier composition of cases differs from that of contacts

Because higher-risk contacts are more likely to become cases, the tier distribution among realised infections is risk-weighted relative to the tier distribution among contacts:

```
    g_l = f_l · r_l / r̄
```

Any quantity that depends on the tier in which a case was infected must be averaged over `g`, not over `f`. In a representative parameterisation with a five-fold gradient in relative risk and tier frequency falling as risk rises, a tracing programme covering the two highest-risk tiers reaches 30% of contacts but 48% of cases. The distinction is carried through the reproduction-number calculation described below.

### Contact tracing

Every contact — not only those that transmit — is independently marked as traced with probability

```
    P(traced | tier l, time t) = c(t) · τ_l
```

where `c(t)` is a time-varying programme coverage lever and `τ_l` is the tier's fixed traceability. This follows the coverage × efficacy factorisation used throughout the model's non-pharmaceutical interventions, in which a time-varying coverage curve is pre-specified from the literature and a scalar efficacy is fitted. Drawing tracing status for every contact, rather than only for contacts that transmit, means the contact log carries the true programme denominator: the number of people a tracing team would have had to enumerate, not merely the number who went on to be infected.

Traced status travels with any contact that becomes a case, and is recorded on the transmission tree. Seed cases have no originating contact and are therefore never traced.

### Effect of tracing on the probability and timing of hospitalisation

The model contains no pre-admission isolation state. The only transmission-reducing condition an infectious individual can occupy is post-admission hospital quarantine, whose efficacy is the mixture

```
    q(t) = ρ(t) · e_ETU + (1 − ρ(t)) · e_hosp
```

over the time-varying proportion `ρ(t)` of admissions managed in an Ebola treatment unit rather than a general hospital, with fixed efficacies `e_ETU` and `e_hosp`. Contact tracing therefore acts on transmission exclusively by delivering cases into that state sooner and more often. It does so through two explicit mechanisms and one implicit one.

**Timing.** A traced case is admitted a flat `d_T` days after symptom onset. This value **caps** rather than replaces the case's own drawn onset-to-admission delay,

```
    D_i = min(D_i^drawn, d_T)   for traced cases
```

so that tracing can only bring an admission forward and never postpone one. A value above the untraced delay distribution would otherwise render an entire tracing scenario a silent no-op; the model therefore warns at the start of a run if `d_T` is not clearly below that distribution.

**Probability.** A traced case may be assigned an absolute probability of hospitalisation, replacing the untraced value outright, or a multiplier on it, capped at one. Both parameterisations are provided because they answer different questions and neither subsumes the other: against a fitted, time-varying baseline hospitalisation curve, a multiplier cannot express a scenario statement of the form "traced cases are hospitalised with probability 0.9", while an absolute probability discards the shape of the fitted curve. Supplying both simultaneously is an error.

**Implicit effect.** Earlier admission raises the *realised* hospitalisation rate independently of either parameter above. An individual is recorded as hospitalised only if their admission time precedes their community outcome time; shortening the onset-to-admission delay therefore increases the proportion of potential admissions that are actually realised, even with the admission probability held fixed.

### Calibration of the basic reproduction number

The contact structure is calibrated through the same single-type, general-population-dominant approximation used previously, evaluated at *t* = 0:

```
    R0_direct  = m_g · p_g · r̄_g · D
    R0_funeral = m_f · p_f · r̄_f · F
    D = 1 − q(0) · Q_g
    F = π_comm · [1 − ε · (1 − u_comm)] + π_hosp · [1 − ε · (1 − u_hosp)]
```

where `Q_g` is the expected fraction of an individual's generation-time mass falling after admission, `π` are the probabilities of dying in the community and in hospital, `u` the corresponding probabilities of an unsafe funeral, and `ε` the efficacy of a safe burial.

Three properties of this construction are worth stating. Contact overdispersion does not appear, because thinning a Negative Binomial leaves its mean unchanged; superspreading can therefore be varied without disturbing a calibration. `Q_g` is raised both by contact tracing, which moves admission earlier, and by the removal of presymptomatic transmission, which deletes the earliest generation-time mass; the two are strongly synergistic, in one test reducing R0 separately by 0.26 and 0.31 but jointly by 0.92. Because whether a case is traced depends on the tier in which it was infected, the Monte Carlo estimation of `Q_g` draws tiers from the case-weighted distribution `g`, not from the contact fractions `f`.

The calculation is separated into an efficacy-**independent** component, computed once by Monte Carlo and cached, and cheap closed-form multipliers recomputed per parameter set, so that an approximate Bayesian computation loop incurs the Monte Carlo cost only once. The relation is inverted to solve for the baseline risks that deliver a target R0 at a specified funeral share of transmission; targets that cannot be attained raise an error naming the achievable ceiling rather than silently clipping a tier's transmission probability.

The funeral share is an **input** to this inversion, consistent with prior practice, and should not be back-calculated from the model's default offspring means. Those defaults are placeholders that any calibrated analysis overwrites, and because `F` embeds both the case fatality risk and safe-burial thinning, the share they imply (approximately 0.09) is substantially below the value typically assumed (0.25) and is not a calibrated quantity.

### Presymptomatic transmission

Contact times and incubation periods are drawn from independent distributions, so a contact may occur before the infector develops symptoms. This behaviour predates the present revision but was not previously visible or controllable. The resulting presymptomatic share is not a parameter but an emergent consequence of the two distributions; for an individual it is

```
    F_GT(T_i^incub) / F_GT(T_i^outcome)
```

where `F_GT` is the generation-time cumulative distribution function. Averaged over simulated natural histories at the parameter values used here — incubation period Gamma with mean 8.5 days and standard deviation 4.5 days, generation time Gamma with mean 15.4 days and shape 2.5 — this yields **36.5%** of transmission occurring before symptom onset, and it ranges from roughly a quarter to a half across plausible filovirus parameter combinations.

This places a ceiling on what any intervention keyed to symptom onset can achieve, contact tracing included, and it is the reason tracing of the highest-risk tier alone produces little reduction in outbreak size in the illustrative analyses. We therefore added a diagnostic that estimates the share once per parameter set and warns when it is large, and a switch that removes presymptomatic transmission entirely by setting the lower truncation bound of the contact-time distribution to the infector's incubation period. The diagnostic saves and restores the random number generator state, so enabling it cannot perturb a simulated trajectory. Note that removing presymptomatic transmission lengthens the realised generation time, since the distribution is being conditioned rather than reshaped; fitted generation-time parameters consequently carry a slightly different meaning afterwards.

### Model outputs and run diagnostics

In addition to the transmission tree, the model now returns a **contact log** with one row per contact generated over the whole run, including contacts that never became infections. Each row carries the infector, the recipient's identifier where transmission occurred, the contact's setting and timing, its risk tier and relative risk, its realised transmission probability, whether it was traced, and — for contacts that did not become infections — the reason, distinguishing failure to transmit from each intervention layer in turn. This is the denominator required for any analysis of tracing or prophylaxis coverage. Because the log is approximately one row per contact rather than per infection, it can be disabled for large calibration runs; doing so consumes no randomness and returns a bit-identical transmission tree.

The model also now reports **why each run terminated**, distinguishing natural extinction from exhaustion of the final-size cap and from exhaustion of the susceptible pool. The latter two return censored final sizes that measure the cap rather than transmission, and are indistinguishable from a controlled outbreak if only the number of cases is examined. A run terminated by the cap raises a warning to this effect. This diagnostic exists because comparisons of mean final size across scenarios that differ in their censoring rate compare censoring, not epidemiology.

---

## 2. Supplementary Tables

**Supplementary Table S1. Parameters introduced or redefined by the contact structure.** All transmission-side parameters are specified independently for the general-population, health-care-worker, and funeral routes.

| Parameter | Interpretation | Notes |
|---|---|---|
| `mn_contacts_*` | Mean number of contacts generated per infectious individual | Replaces the former offspring mean. Scalar or function of calendar time |
| `overdisp_contacts_*` | Negative Binomial overdispersion of the contact count | Replaces the former offspring overdispersion. Does not enter R0 |
| `baseline_risk_*` | Per-contact transmission probability of a reference-tier contact | Under the default reference convention, the transmission probability of the highest-risk contact |
| `contact_risk_*` | Risk-tier structure for the route | Falls back to a single shared structure where not specified per route |
| — `fractions` | Share of contacts in each tier | Must sum to 1 |
| — `relative_risk` | Transmission risk of each tier relative to the reference | Normalised so the reference tier equals 1 |
| — `trace_prob` | Probability a tier-*l* contact is traceable | Multiplied by the time-varying programme coverage |
| — `reference` | Tier anchoring the relative-risk scale | Defaults to the highest-risk tier |
| `trace_coverage` | Programme-level tracing coverage | Scalar or function of calendar time |
| `onset_to_hospitalisation_traced` | Flat onset-to-admission delay for traced cases, in days | Caps rather than replaces each case's drawn delay |
| `prob_hospitalised_traced` | Absolute probability of hospitalisation for traced cases | Mutually exclusive with the multiplier below |
| `prob_hospitalised_multiplier_traced` | Multiplier on the untraced hospitalisation probability | Capped at 1 |
| `presymptomatic_transmission` | Whether contacts may occur before the infector's symptom onset | Removal is by exact truncation of the contact-time distribution |
| `r0_target`, `r0_prop_funeral` | Target R0 and funeral share, inverted to baseline risks at run start | Solved values reported with the run |

**Supplementary Table S2. Illustrative risk-tier structure and its derived quantities.** Five tiers with relative risks spanning a five-fold range, with tier frequency falling as risk rises, such that the lowest-risk tier is exactly five times as common as the highest. Note that the composition of cases is markedly flatter than the composition of contacts: rare, high-risk exposures are substantially over-represented among infections.

| Tier | Share of contacts, *f_l* | Relative risk, *r_l* | Share of cases, *g_l* |
|---|---|---|---|
| 1 | 0.333 | 0.2 | 0.143 |
| 2 | 0.267 | 0.4 | 0.229 |
| 3 | 0.200 | 0.6 | 0.257 |
| 4 | 0.133 | 0.8 | 0.229 |
| 5 (reference) | 0.067 | 1.0 | 0.143 |

Mean relative risk `r̄` = 0.467. With mean contacts of 15 for the direct routes and 20 for the funeral route, and a target R0 of 1.35 at a funeral share of 0.25, the inversion returns baseline per-contact transmission probabilities of 0.150 for the direct routes and 0.073 for the funeral route.

---

## 3. Known limitations of the current implementation

Four limitations bear on the interpretation of tracing analyses in particular.

**Contact number does not depend on the duration of infectiousness.** The contact count is drawn before the individual's outcome time is consulted, and contact times are then truncated into the surviving window. An individual who dies on day 10 therefore generates the same expected number of contacts as one who survives to day 60; the contacts are compressed into the shorter interval rather than reduced in number. Since admission shortens the time to outcome by approximately a fifth in this parameterisation, the transmission-reducing benefit of earlier admission is represented only through relocation of contacts into the quarantined window, and not at all through shortening of the infectious period. Estimated tracing benefits are consequently conservative. A rate-based formulation, in which the contact mean scales with expected time to outcome, would remove this limitation at the cost of requiring recalibration and a re-derivation of the reproduction-number approximation.

**Contact tracing reshapes the distribution of burden as well as its total.** Earlier admission relocates transmission from the community into the hospital, where contacts are substantially more likely to be health-care workers. Total transmission falls under tracing in every scenario examined, but the health-care-worker share of infections can rise, and does so most where Ebola treatment unit capacity and personal protective equipment coverage are lowest — that is, earliest in a response. Analyses of tracing should report health-care-worker burden separately from total burden. This interacts with the limitation above: because contacts are relocated rather than removed, the magnitude of this reallocation is likely to be overstated.

**The tier profile of tracing is fixed in shape.** Programme coverage varies with time but the relative traceability of the tiers does not, so the model can represent a programme scaling up or down but not one that broadens its reach from household contacts outward as capacity accumulates — which is the realistic trajectory, and the one that would plausibly avoid the early reallocation of burden described above.

**Susceptible depletion terminates a run rather than attenuating transmission.** The susceptible pool is decremented by each infection but does not feed back into the per-contact transmission probability, so an epidemic runs at undiminished force until the pool is exhausted and then stops. This has no effect at the population sizes used here, where the pool cannot approach exhaustion, but it would misrepresent any outbreak modelled in a small closed population. The contact structure makes the remedy straightforward — an additional multiplicative susceptible-fraction term in the per-contact transmission probability — at the cost of converting the reproduction-number approximation from an R0 into an explicitly *t* = 0 quantity.
