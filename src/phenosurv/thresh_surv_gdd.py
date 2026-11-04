"""Threshold survival model with growing-degree-day (GDD) forcing."""

from importlib.resources import files

# pyrefly: ignore[missing-import]
import cmdstanpy


class ThreshSurvGDD:
    def __init__(self):
        self.stan_file = files("phenosurv") / "stan" / "thresh-surv-gdd.stan"
        self.model = cmdstanpy.CmdStanModel(
            stan_file=self.stan_file,
            )

    def print_code(self):
        """Print the model's Stan code."""
        print(self.stan_file.read_text())
