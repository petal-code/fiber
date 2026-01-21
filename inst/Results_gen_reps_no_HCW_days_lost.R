## ================================================================
## 200 stochastic repeats: branching_process_main() -> summarise_output()
## Uses the same parameter values as previous tests, but runs many seeds.
### note as of yet hcw_loss_function is not processed, this needs to be updated so that HCW loss metrics are included in
### summarise_output - once they are - this code should accommodate for that.
### note, HCW_per_capita error still needs to be investigated.
## ================================================================

source("R/branching_process_main.R")
source("R/summarise_output.R")
source("helper_functions.R")
source("R/offspring_function_genPop.R")
source("R/offspring_function_hcw.R")
source("R/offspring_function_funeral.R")
source("R/complete_offspring_info.R")


## ---- 1) Natural history functions (identical to previous tests)
incubation_period <- function(n) rgamma(n, shape = 2.5, rate = 0.5)
onset_to_hospitalisation <- function(n) rgamma(n, shape = 2, rate = 0.5)
onset_to_death <- function(n) rgamma(n, shape = 2, rate = 0.25)
onset_to_recovery <- function(n) rgamma(n, shape = 2, rate = 0.25)
hospitalisation_to_death <- function(n) rgamma(n, shape = 2, rate = 0.4)
hospitalisation_to_recovery <- function(n) rgamma(n, shape = 2, rate = 0.4)

## ---- 2) Parameters
mn_offspring_genPop        <- 4
overdisp_offspring_genPop  <- 0.5
Tg_shape_genPop            <- 2
Tg_rate_genPop             <- 0.5

mn_offspring_hcw           <- 1.8
overdisp_offspring_hcw     <- 0.5
Tg_shape_hcw               <- 2
Tg_rate_hcw                <- 0.5

mn_offspring_funeral       <- 3
overdisp_offspring_funeral <- 0.5
Tg_shape_funeral           <- 3
Tg_rate_funeral            <- 1

prob_symptomatic              <- 0.8
prob_hospitalised_hcw         <- 0.7
prob_hospitalised_genPop      <- 0.6
prob_death_comm               <- 0.5
prob_death_hosp               <- 0.4

prob_hcw_cond_genPop_comm     <- 0.05
prob_hcw_cond_genPop_hospital <- 0.30

prob_hcw_cond_hcw_comm        <- 0.10
prob_hcw_cond_hcw_hospital    <- 0.50

prob_hospital_cond_hcw_preAdm <- 0.7
ppe_efficacy_hcw              <- 0.5
hospital_quarantine_efficacy  <- 0.8

p_unsafe_funeral_comm_hcw      <- 0.4
p_unsafe_funeral_hosp_hcw      <- 0.4
p_unsafe_funeral_comm_genPop   <- 0.5
p_unsafe_funeral_hosp_genPop   <- 0.5
safe_funeral_efficacy          <- 0.8

prob_hcw_cond_funeral_hcw      <- 0.10
prob_hcw_cond_funeral_genPop   <- 0.05

population       <- 10000
hcw_per_capita   <- 0.1     # override to a "working" value for now
check_final_size <- 5000
initial_immune   <- 0
seeding_cases    <- 10

## ---- 3) Fixed argument list (everything except seed) ----
base_args <- list(
  ## Transmission
  mn_offspring_genPop        = mn_offspring_genPop,
  overdisp_offspring_genPop  = overdisp_offspring_genPop,
  Tg_shape_genPop            = Tg_shape_genPop,
  Tg_rate_genPop             = Tg_rate_genPop,
  mn_offspring_hcw           = mn_offspring_hcw,
  overdisp_offspring_hcw     = overdisp_offspring_hcw,
  Tg_shape_hcw               = Tg_shape_hcw,
  Tg_rate_hcw                = Tg_rate_hcw,
  mn_offspring_funeral       = mn_offspring_funeral,
  overdisp_offspring_funeral = overdisp_offspring_funeral,
  Tg_shape_funeral           = Tg_shape_funeral,
  Tg_rate_funeral            = Tg_rate_funeral,

  ## Natural history
  incubation_period           = incubation_period,
  onset_to_hospitalisation    = onset_to_hospitalisation,
  onset_to_death              = onset_to_death,
  onset_to_recovery           = onset_to_recovery,
  hospitalisation_to_death    = hospitalisation_to_death,
  hospitalisation_to_recovery = hospitalisation_to_recovery,

  ## Disease severity and healthcare seeking
  prob_symptomatic         = prob_symptomatic,
  prob_hospitalised_hcw    = prob_hospitalised_hcw,
  prob_hospitalised_genPop = prob_hospitalised_genPop,
  prob_death_comm          = prob_death_comm,
  prob_death_hosp          = prob_death_hosp,

  ## genPop → {genPop, HCW}
  prob_hcw_cond_genPop_comm     = prob_hcw_cond_genPop_comm,
  prob_hcw_cond_genPop_hospital = prob_hcw_cond_genPop_hospital,

  ## HCW → {genPop, HCW}
  prob_hcw_cond_hcw_comm     = prob_hcw_cond_hcw_comm,
  prob_hcw_cond_hcw_hospital = prob_hcw_cond_hcw_hospital,

  ## HCW setting model
  prob_hospital_cond_hcw_preAdm = prob_hospital_cond_hcw_preAdm,
  ppe_efficacy_hcw              = ppe_efficacy_hcw,
  hospital_quarantine_efficacy  = hospital_quarantine_efficacy,

  ## Funeral occurrence
  p_unsafe_funeral_comm_hcw     = p_unsafe_funeral_comm_hcw,
  p_unsafe_funeral_hosp_hcw     = p_unsafe_funeral_hosp_hcw,
  p_unsafe_funeral_comm_genPop  = p_unsafe_funeral_comm_genPop,
  p_unsafe_funeral_hosp_genPop  = p_unsafe_funeral_hosp_genPop,
  safe_funeral_efficacy         = safe_funeral_efficacy,

  ## HCW vs genPop at funeral
  prob_hcw_cond_funeral_hcw    = prob_hcw_cond_funeral_hcw,
  prob_hcw_cond_funeral_genPop = prob_hcw_cond_funeral_genPop,

  ## Misc
  tf                  = Inf,
  population          = population,
  hcw_per_capita      = hcw_per_capita,
  check_final_size    = check_final_size,
  initial_immune      = initial_immune,
  seeding_cases       = seeding_cases,
  susceptible_deplete = FALSE
)

## ---- 4) One-run sanity check
## summarise_output expects (tdf, subset, sim_info).
sim1 <- do.call(branching_process_main, c(base_args, list(seed = 123)))
summary1 <- summarise_output(tdf = sim1$tdf, subset = "realised_subset", sim_info = sim1$sim_info)
print(summary1[c("n_cases_total","n_cases_genPop","n_cases_HCW","n_funeral","prop_funeral",
                 "attack_rate_overall","attack_rate_genPop","hcw_attack_rate","deaths_per_1000_pop")])

## ---- 5) 200 repeats: ONLY the seed changes ----
run_one <- function(seed) {
  out  <- do.call(branching_process_main, c(base_args, list(seed = seed)))
  summ <- summarise_output(tdf = out$tdf, subset = "realised_subset", sim_info = out$sim_info)
  ## return one row data.frame
  data.frame(seed = seed, failed = FALSE, as.data.frame(summ), check.names = FALSE)
}

seeds <- 1:200

res_list <- lapply(seeds, function(s) {
  tryCatch(
    run_one(s),
    error = function(e) data.frame(seed = s, failed = TRUE, error = conditionMessage(e))
  )
})

res_200 <- do.call(rbind, res_list)

## ---- 6) Save outputs ----
write.csv(res_200, "gilead_repeats_200_raw.csv", row.names = FALSE)

## ---- 7) Quick aggregate summary (optional) ----
ok <- res_200$failed == FALSE

p_gt_20  <- mean(res_200$n_cases_total[ok] > 20,  na.rm = TRUE)
p_gt_100 <- mean(res_200$n_cases_total[ok] > 100, na.rm = TRUE)

summary_tbl <- data.frame(
  reps_ok = sum(ok),
  p_gt_20 = p_gt_20,
  p_gt_100 = p_gt_100,
  med_cases = median(res_200$n_cases_total[ok], na.rm = TRUE),
  lo_cases  = as.numeric(quantile(res_200$n_cases_total[ok], 0.025, na.rm = TRUE)),
  hi_cases  = as.numeric(quantile(res_200$n_cases_total[ok], 0.975, na.rm = TRUE)),
  med_deaths = median(res_200$n_deaths_total[ok], na.rm = TRUE),
  med_hcw_cases = median(res_200$n_cases_HCW[ok], na.rm = TRUE),
  med_duration_days = median(res_200$outbreak_duration_days[ok], na.rm = TRUE)
)

write.csv(summary_tbl, "gilead_repeats_200_summary.csv", row.names = FALSE)
print(summary_tbl)



