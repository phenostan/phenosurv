"""Threshold survival model with growing-degree-day (GDD) forcing."""

from importlib.resources import files

# pyrefly: ignore[missing-import]
import cmdstanpy
import pandas as pd


class ThreshSurvGDD:
    def __init__(self):
        self.stan_file = files("phenosurv") / "stan" / "thresh-surv-gdd.stan"
        self.model = cmdstanpy.CmdStanModel(
            stan_file=self.stan_file,
            )
        self.fit_obj = None

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

    def summary(self) -> pd.DataFrame:
        """
        Generates a summary table of the posterior distribution

        Returns: 
        pd.DataFrame: Summary statistics (Mean, MCSE, StdDev, qantiles, 
            N_eff, R_hat).
        """
        self._check_is_fitted()
        return self.fit_obj.summary()

        


