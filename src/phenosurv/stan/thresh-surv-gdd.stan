/* MODIFIED (THRESHOLDING) SURVIVAL MODEL w/ LINEAR (GDD) FORCING -- THREADED   */
/* simple growing-degree-day hinge:                                            */
/*                                                                             */
/*     forcing(T) = max(T - T_base, 0)        (FU / day = degree-days)         */
/*                                                                             */
/* ########################################################################### */ 
// #
// # The code is copyright 2024 Michael Betancourt and licensed under the
// # new BSD (3-clause) license:
// #  https://opensource.org/licenses/BSD-3-Clause
// #
// # For more information see https://github.com/betanalpha/mcmc_diagnostics.
// #
// ############################################################################

functions {

  // Linear (growing-degree-day) forcing increment.
  real forcing(real T, real T_base) {
    return fmax(T - T_base, 0.0);
  }

  // log of the instantaneous forcing rate at the event day (the dPsi/dt term).
  real log_forcing(real T, real T_base) {
    if (T <= T_base) return negative_infinity();
    return log(T - T_base);
  }

  // PARALLEL PARTIAL SUM FUNCTION
  // Computes the log-likelihood contribution of a subset of site-years.
  real partial_sum(
    array[] int site_year_slice,
    int start,
    int end,
    array[] real obs_temps,
    array[] int temp_start_idxs,
    array[] int N_days,
    array[] int precursor_events,
    array[] int events,
    array[] int genotype,
    array[] int event_start_idxs,
    array[] int N_events,
    real T_base,
    vector Psi0,
    real sigma
  ) {
    real lp = 0;

    for (i in 1:size(site_year_slice)) {
      int site_year_idx = site_year_slice[i];

      int temp_start_idx = temp_start_idxs[site_year_idx];
      int temp_end_idx   = temp_start_idxs[site_year_idx] + N_days[site_year_idx] - 1;
      array[N_days[site_year_idx]] real local_temps = obs_temps[temp_start_idx:temp_end_idx];

      int event_start_idx = event_start_idxs[site_year_idx];
      int event_end_idx   = event_start_idxs[site_year_idx] + N_events[site_year_idx] - 1;
      array[N_events[site_year_idx]] int local_precursor_events = precursor_events[event_start_idx:event_end_idx];
      array[N_events[site_year_idx]] int local_events           = events[event_start_idx:event_end_idx];

      for (n in 1:N_events[site_year_idx]) {
        int trial = event_start_idx + n - 1;
        int g = genotype[trial];
        int duration = local_events[n] - local_precursor_events[n];
        vector[duration] daily_forcings;

        for (j in 1:duration) {
          int day = local_precursor_events[n] + j - 1;
          daily_forcings[j] = forcing(local_temps[day], T_base);
        }

        real Psi        = sum(daily_forcings);
        real log_dPsidt = log_forcing(local_temps[local_events[n]], T_base);
        lp += logistic_lpdf(Psi | Psi0[g], sigma) + log_dPsidt;
      }
    }

    return lp;
  }
}

data {
  int N_temp_grid;                          // Size of grid of relevant temperatures
  array[N_temp_grid] real temp_grid;        // Grid of relevant temperatures (C)

  int N_site_years;                         // Number of year-site instances
  int N_obs_temps;                          // Number of year-site days
  array[N_obs_temps] real obs_temps;        // Recorded temperature (C) for each day

  array[N_site_years] int temp_start_idxs;
  array[N_site_years] int N_days;

  int N;                                    // Number of phenological events
  array[N] int precursor_events;            // Day transition starts (<= 366)
  array[N] int events;                      // Day transition ends   (<= 366)

  int N_geno;                               // Number of unique genotypes
  array[N] int genotype;                    // Genotype ID for each event

  array[N_site_years] int event_start_idxs;
  array[N_site_years] int N_events;

  // Control parameters for threading
  int<lower=0,upper=1> use_threading;       // 1 = use threading, 0 = don't
  int<lower=1> grainsize;                   // Minimum site-years per thread (1 = auto)
}

transformed data {
  // Create array of site-year indices for reduce_sum
  array[N_site_years] int site_year_indices;
  for (i in 1:N_site_years) {
    site_year_indices[i] = i;
  }
}

parameters {
  real T_base;                              // Base temperature for forcing (C)

  // Genotype-varying forcing requirement F* (= Psi0), non-centered.
  real Psi0_bar;
  vector[N_geno] z_psi0;
  real<lower=0> tau_psi0;

  real<lower=0> sigma;                      // Threshold softness scale (FU)
}

transformed parameters {
  vector[N_geno] Psi0 = Psi0_bar + z_psi0 * tau_psi0;
}

model {
  // Priors --------------------------------------------------------------
  // Base temperature: sunflower GDD base commonly ~6.7 C (Kalyar et al. 2014).
  T_base ~ normal(4.8, 3 / 2.32);           // ~ [1.8, 7.8] C at 98%

  // F* in degree-day units (NOT comparable in scale to the W&E Psi0).
  // ~70 d flowering window at ~13 C above base => F* on the order of ~900 DD.
  Psi0_bar ~ normal(900, 250 / 2.32);
  z_psi0   ~ normal(0, 1);
  tau_psi0 ~ normal(0, 150);

  sigma ~ normal(0, 150 / 2.57);

  // Likelihood with optional threading ---------------------------------
  if (use_threading == 1) {
    target += reduce_sum(
      partial_sum,
      site_year_indices,
      grainsize,
      obs_temps,
      temp_start_idxs,
      N_days,
      precursor_events,
      events,
      genotype,
      event_start_idxs,
      N_events,
      T_base,
      Psi0,
      sigma
    );
  } else {
    for (i in 1:N_site_years) {
      int temp_start_idx = temp_start_idxs[i];
      array[N_days[i]] real local_temps = obs_temps[temp_start_idx:temp_start_idx + N_days[i] - 1];

      int event_start_idx = event_start_idxs[i];
      array[N_events[i]] int local_precursor_events = precursor_events[event_start_idx:event_start_idx + N_events[i] - 1];
      array[N_events[i]] int local_events           = events[event_start_idx:event_start_idx + N_events[i] - 1];

      for (n in 1:N_events[i]) {
        int trial = event_start_idx + n - 1;
        int g = genotype[trial];
        int duration = local_events[n] - local_precursor_events[n];
        vector[duration] daily_forcings;

        for (j in 1:duration) {
          int day = local_precursor_events[n] + j - 1;
          daily_forcings[j] = forcing(local_temps[day], T_base);
        }

        real Psi        = sum(daily_forcings);
        real log_dPsidt = log_forcing(local_temps[local_events[n]], T_base);
        target += logistic_lpdf(Psi | Psi0[g], sigma) + log_dPsidt;
      }
    }
  }
}

generated quantities {
  vector[N_temp_grid] forcings;             // forcing curve on the grid (geno-independent)
  array[N] int pred_events;                 // retrodicted flowering DOY
  vector[N] log_lik;                        // pointwise log-lik for PSIS-LOO

  // Forcing curve (shape identical across genotypes here; requirement varies).
  for (t in 1:N_temp_grid)
    forcings[t] = forcing(temp_grid[t], T_base);

  {
    for (i in 1:N_site_years) {
      int temp_start_idx = temp_start_idxs[i];
      array[N_days[i]] real local_temps = obs_temps[temp_start_idx:temp_start_idx + N_days[i] - 1];

      int event_start_idx = event_start_idxs[i];
      array[N_events[i]] int local_precursor_events = precursor_events[event_start_idx:event_start_idx + N_events[i] - 1];
      array[N_events[i]] int local_events           = events[event_start_idx:event_start_idx + N_events[i] - 1];

      int first_precursor_event = min(local_precursor_events);
      vector[N_days[i]] daily_forcings = rep_vector(0, N_days[i]);
      vector[N_days[i]] accum_forcings;
      for (t in first_precursor_event:N_days[i])
        daily_forcings[t] = forcing(local_temps[t], T_base);
      accum_forcings = cumulative_sum(daily_forcings);

      for (n in 1:N_events[i]) {
        int trial = event_start_idx + n - 1;
        int g = genotype[trial];
        int precursor = local_precursor_events[n];
        int event_day = local_events[n];

        // pointwise log_lik -- EXACTLY mirrors the model{} increment:
        //   Psi = sum_{day=precursor}^{event_day-1} forcing(T_day)
        //       = accum[event_day-1] - accum[precursor-1]
        real Psi_obs    = accum_forcings[event_day - 1]
                        - (precursor > 1 ? accum_forcings[precursor - 1] : 0);
        real log_dPsidt = log_forcing(local_temps[event_day], T_base);
        log_lik[trial]  = logistic_lpdf(Psi_obs | Psi0[g], sigma) + log_dPsidt;

        // retrodicted event day
        real init_forcing = accum_forcings[precursor];
        real l = logistic_rng(Psi0[g], sigma);
        pred_events[trial] = N_days[i];
        for (d in precursor:N_days[i]) {
          if (accum_forcings[d] - init_forcing >= l) {
            pred_events[trial] = d;
            break;
          }
        }
      }
    }
  }
}
