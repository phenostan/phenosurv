"""Simulate a small dataset from the thresh-surv-gdd model for the test suite.

Writes two files to tests/data/ (or --out-dir):

    testdata.json        Stan data for thresh-surv-gdd.stan (ready for cmdstanpy)
    testdata_truth.json  the parameter values the data were simulated from

Run through pixi:

    pixi run simulate-testdata
    pixi run simulate-testdata --seed 7 --n-site-years 10

Generative process (mirrors the model block of thresh-surv-gdd.stan):

    T_base, Psi0_bar, tau_psi0, sigma          fixed "true" values
    Psi0[g] = Psi0_bar + z[g] * tau_psi0       z[g] ~ normal(0, 1)
    obs_temps                                   seasonal sine + daily noise, per site-year
    threshold l ~ logistic(Psi0[g], sigma)      one draw per event
    event = first day t >= precursor + 1 with
            sum_{d = precursor}^{t} max(T_d - T_base, 0) >= l

That event rule is the discrete form of the likelihood
logistic_lpdf(Psi | Psi0, sigma) + log(dPsi/dt), where
Psi = sum_{d = precursor}^{event - 1} forcing(T_d): the threshold is crossed
during the event day, so forcing on the event day is always > 0.
"""

import argparse
import json
from pathlib import Path

import numpy as np

OUT_DIR = Path(__file__).resolve().parent.parent / "tests" / "data"

TRUTH = {
    "T_base": 5.0,       # C
    "Psi0_bar": 900.0,   # degree-days
    "tau_psi0": 80.0,    # degree-days
    "sigma": 40.0,       # degree-days
}

N_DAYS = 365


def simulate_temps(rng):
    """One site-year of daily mean temperature (C): seasonal sine plus noise."""
    doy = np.arange(1, N_DAYS + 1)
    mean = rng.normal(10.0, 2.0)          # site-year annual mean
    amplitude = rng.normal(12.0, 1.5)     # seasonal swing
    seasonal = mean + amplitude * np.sin(2 * np.pi * (doy - 110) / 365)
    return seasonal + rng.normal(0.0, 3.0, size=N_DAYS)


def event_day(temps, precursor, threshold, T_base):
    """First 1-based day the accumulated forcing from `precursor` reaches `threshold`."""
    forcing = np.maximum(temps[precursor - 1:] - T_base, 0.0)
    crossed = np.nonzero(np.cumsum(forcing) >= threshold)[0]
    return None if crossed.size == 0 else precursor + int(crossed[0])


def simulate(rng, n_site_years, n_geno, reps):
    T_base = TRUTH["T_base"]
    z_psi0 = rng.normal(0.0, 1.0, size=n_geno)
    Psi0 = TRUTH["Psi0_bar"] + z_psi0 * TRUTH["tau_psi0"]

    obs_temps, temp_start_idxs, N_days = [], [], []
    precursor_events, events, genotype = [], [], []
    event_start_idxs, N_events = [], []

    for _ in range(n_site_years):
        temps = simulate_temps(rng)
        temp_start_idxs.append(len(obs_temps) + 1)
        N_days.append(N_DAYS)
        obs_temps.extend(temps.tolist())

        event_start_idxs.append(len(events) + 1)
        planting = int(rng.integers(110, 151))          # site-year planting DOY
        n = 0
        for g in range(n_geno):
            for _ in range(reps):
                precursor = planting
                while True:
                    threshold = rng.logistic(Psi0[g], TRUTH["sigma"])
                    if threshold <= 0:
                        continue
                    day = event_day(temps, precursor, threshold, T_base)
                    if day is not None and day > precursor:
                        break
                precursor_events.append(precursor)
                events.append(day)
                genotype.append(g + 1)                  # Stan is 1-indexed
                n += 1
        N_events.append(n)

    temp_grid = np.linspace(-5.0, 35.0, 41)
    data = {
        "N_temp_grid": len(temp_grid),
        "temp_grid": temp_grid.round(6).tolist(),
        "N_site_years": n_site_years,
        "N_obs_temps": len(obs_temps),
        "obs_temps": np.round(obs_temps, 3).tolist(),
        "temp_start_idxs": temp_start_idxs,
        "N_days": N_days,
        "N": len(events),
        "precursor_events": precursor_events,
        "events": events,
        "N_geno": n_geno,
        "genotype": genotype,
        "event_start_idxs": event_start_idxs,
        "N_events": N_events,
        "use_threading": 0,
        "grainsize": 1,
    }
    truth = {**TRUTH, "z_psi0": z_psi0.tolist(), "Psi0": Psi0.tolist()}
    return data, truth


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, default=20261103)
    parser.add_argument("--n-site-years", type=int, default=4)
    parser.add_argument("--n-geno", type=int, default=3)
    parser.add_argument("--reps", type=int, default=2, help="events per genotype per site-year")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    data, truth = simulate(rng, args.n_site_years, args.n_geno, args.reps)
    truth["seed"] = args.seed

    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "testdata.json").write_text(json.dumps(data) + "\n")
    (args.out_dir / "testdata_truth.json").write_text(json.dumps(truth, indent=2) + "\n")
    print(f"Wrote {data['N']} events across {data['N_site_years']} site-years to {args.out_dir}")


if __name__ == "__main__":
    main()
