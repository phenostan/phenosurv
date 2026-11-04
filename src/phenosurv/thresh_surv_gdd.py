"""Threshold survival model with growing-degree-day (GDD) forcing."""

import json
import os
from importlib.resources import files

# pyrefly: ignore[missing-import]
import cmdstanpy
import numpy as np
import pandas as pd


class ThreshSurvGDD:
    def __init__(self):
        self.stan_file = files("phenosurv") / "stan" / "thresh-surv-gdd.stan"
        self.model = cmdstanpy.CmdStanModel(
            stan_file=self.stan_file,
            )
        self.fit_obj = None
        self.data = None

    def print_code(self):
        """Print the model's Stan code."""
        print(self.stan_file.read_text())

    def fit(
        self, 
        data: dict,
        chains: int = 4,
        iter_warmup: int = 1000,
        iter_sampling: int = 1000,
        **kwargs):
        """
        Fits the Stan model to the provided `data` using the No-U-Turn-Sampler 
        variant of Hamiltonian Monte Carlo via `cmdstanpy`. 

        Args:
            data (dict): Dictionary of data inputs required by the Stan model.
            chains (int): Number of MCMC chains.
            iter_warmup (int): Number of warmup iterations per chain.
            iter_sampling (int): Number of sampling iterations per chain.
            **kwargs: Additional arguments passed to cmdstanpy.CmdStanModel.sample()

        Returns:
            self: Returns the instance so methods can be chained. 
        """

        self.data = data
        self.fit_obj = self.model.sample(
            data=data,
            chains=chains,
            iter_warmup=iter_warmup,
            iter_sampling=iter_sampling,
            **kwargs
        )
        return self
    
    def _check_is_fitted(self):
        """Internal helper to ensure the model has been fitted...
        """
        if self.fit_obj == None:
            raise RuntimeError("Model has not been fitted. Call .fit() first...")

    def diagnose(self) -> str:
        """
        Runs cmdtan's built-in diagnostics for HMC/NUTS, etc.

        Returns:
            str: Daignostic output message.
        """
        self._check_is_fitted()
        return self.fit_obj.diagnose()

    def predict(self, data: dict | str | os.PathLike | None = None, **kwargs) -> np.ndarray:
        """
        Predicts event days (day of year) from the fitted posterior.

        For every posterior draw, the model accumulates forcing from each
        event's precursor day and draws a forcing threshold, giving one
        predicted event day per draw. Uses cmdstanpy's standalone generated
        quantities, so the model is not refit.

        Args:
            data (dict | str | PathLike | None): Data in the same format as
                for `fit()` (a dict or a path to a JSON file). `events` may be
                omitted, since it is what is being predicted. Defaults to the
                data the model was fitted to.
            **kwargs: Additional arguments passed to
                cmdstanpy.CmdStanModel.generate_quantities() (e.g. `seed`).

        Returns:
            np.ndarray: Integer array of shape (n_draws, N) of predicted event
                days, one column per event in `data`.
        """
        self._check_is_fitted()
        data = self._load_data(self.data if data is None else data)
        if "events" not in data:
            data["events"] = self._placeholder_events(data)

        gq = self.model.generate_quantities(
            data=data,
            previous_fit=self.fit_obj,
            **kwargs,
        )
        return np.asarray(gq.stan_variable("pred_events")).astype(int)
    def summary(self) -> pd.DataFrame:
        """
        Generates a summary table of the posterior distribution

        Returns: 
        pd.DataFrame: Summary statistics (Mean, MCSE, StdDev, qantiles, 
            N_eff, R_hat).
        """
        self._check_is_fitted()
        return self.fit_obj.summary()

        

    @staticmethod
    def _load_data(data) -> dict:
        """Returns a shallow copy of `data` as a dict, reading it if a JSON path."""
        if isinstance(data, (str, os.PathLike)):
            with open(data) as f:
                return json.load(f)
        return dict(data)

    @staticmethod
    def _placeholder_events(data: dict) -> list[int]:
        """
        Dummy event days for prediction when `events` are unknown.

        The Stan program needs `events` to compute `log_lik`, which is ignored
        by `predict()`; `pred_events` does not depend on them. Each placeholder
        is the day after the precursor, capped at the site-year's last day.
        """
        events = []
        for start, n_events, n_days in zip(
            data["event_start_idxs"], data["N_events"], data["N_days"]
        ):
            for trial in range(start - 1, start - 1 + n_events):
                events.append(min(data["precursor_events"][trial] + 1, n_days))
        return events
