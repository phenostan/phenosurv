"""Tests for the ThreshSurvGDD model class (threshold survival model, GDD forcing)."""

import contextlib
import inspect
import io
import json
import unittest
from importlib.resources import files
from pathlib import Path

# pyrefly: ignore[missing-import]
import cmdstanpy
import numpy as np
import pandas as pd

from phenosurv import ThreshSurvGDD


class TestThreshSurvGDDClass(unittest.TestCase):

    def test_is_a_class(self):
        self.assertTrue(inspect.isclass(ThreshSurvGDD))

    def test_can_instantiate_without_data(self):
        # Inspecting the Stan code should not require data.
        model = ThreshSurvGDD()
        self.assertIsInstance(model, ThreshSurvGDD)

    def test_model_is_compiled_cdmstanmodel(self):
        m = ThreshSurvGDD()
        
        # check self.model compiles and is a CmdStanModel instance
        self.assertIsInstance(m.model, cmdstanpy.CmdStanModel)

        # verify the model has its c++ executeable
        self.assertIsNotNone(m.model.exe_file)


class TestThreshSurvGDDPrintCode(unittest.TestCase):

    def setUp(self):
        self.model = ThreshSurvGDD()

    def _printed_code(self):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            self.model.print_code()
        return buffer.getvalue()

    def test_has_print_code_method(self):
        self.assertTrue(callable(getattr(self.model, "print_code", None)))

    def test_print_code_writes_to_stdout(self):
        self.assertNotEqual(self._printed_code().strip(), "")

    def test_printed_code_contains_stan_program_blocks(self):
        code = self._printed_code()
        for block in (
            "functions {",
            "data {",
            "transformed data {",
            "parameters {",
            "transformed parameters {",
            "model {",
            "generated quantities {",
        ):
            with self.subTest(block=block):
                self.assertIn(block, code)

    def test_printed_code_is_the_gdd_threshold_model(self):
        code = self._printed_code()
        self.assertIn("MODIFIED (THRESHOLDING) SURVIVAL MODEL w/ LINEAR (GDD) FORCING", code)
        for declaration in ("real T_base;", "real Psi0_bar;", "real<lower=0> sigma;"):
            with self.subTest(declaration=declaration):
                self.assertIn(declaration, code)

    def test_printed_code_matches_packaged_stan_file(self):
        packaged = (files("phenosurv") / "stan" / "thresh-surv-gdd.stan").read_text()
        self.assertEqual(self._printed_code().rstrip("\n"), packaged.rstrip("\n"))

class TestThreshSurvGDDFitted(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Runs once for the whole class. Great for slow operations like MCMC."""
        cls.m = ThreshSurvGDD()
        cls.path_test_data = Path("tests", "data", "testdata.json")
        cls.fit_args = {
            "chains": 2,
            "iter_warmup": 300,
            "iter_sampling": 300,
            "seed": 1997,
            "show_progress": False
        }

        cls.m.fit(data=str(cls.path_test_data), **cls.fit_args)

    def test_fitted_is_cmdstanmcmc(self):
        self.assertIsInstance(self.m.fit_obj, cmdstanpy.CmdStanMCMC)

    def test_fitted_diagnoses_returns_string(self):
            output = self.m.diagnose()

            self.assertNotEqual(output.strip(), "")
            
            self.assertIn("transitions treedepth", output.lower())
            self.assertIn("divergences", output.lower())
            self.assertIn("e-bfmi", output.lower())
            self.assertIn("rank-normalized split", output.lower())

    def test_fitted_summary(self):
        summary = self.m.summary()

        self.assertIsInstance(summary, pd.DataFrame)
        self.assertGreater(summary.shape[0], 1)


class TestThreshSurvGDDPredictUnfitted(unittest.TestCase):

    def test_predict_before_fit_raises(self):
        m = ThreshSurvGDD()
        with self.assertRaises(RuntimeError):
            m.predict()


class TestThreshSurvGDDPredict(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.m = ThreshSurvGDD()
        cls.path_test_data = Path("tests", "data", "testdata.json")
        with open(cls.path_test_data) as f:
            cls.data = json.load(f)
        cls.chains = 2
        cls.iter_sampling = 300
        cls.n_draws = cls.chains * cls.iter_sampling
        cls.m.fit(
            data=str(cls.path_test_data),
            chains=cls.chains,
            iter_warmup=300,
            iter_sampling=cls.iter_sampling,
            seed=1997,
            show_progress=False,
        )

    def _first_site_year(self):
        """New data: only the first site-year of the test data, without `events`."""
        d = self.data
        n_days = d["N_days"][0]
        n_events = d["N_events"][0]
        new = dict(d)
        new.update(
            N_site_years=1,
            N_obs_temps=n_days,
            obs_temps=d["obs_temps"][:n_days],
            temp_start_idxs=[1],
            N_days=[n_days],
            N=n_events,
            precursor_events=d["precursor_events"][:n_events],
            genotype=d["genotype"][:n_events],
            event_start_idxs=[1],
            N_events=[n_events],
        )
        del new["events"]
        return new

    def test_has_predict_method(self):
        self.assertTrue(callable(getattr(self.m, "predict", None)))

    def test_predict_defaults_to_training_data(self):
        pred = self.m.predict(seed=1)
        self.assertIsInstance(pred, np.ndarray)
        self.assertEqual(pred.shape, (self.n_draws, self.data["N"]))

    def test_predict_accepts_dict_and_path(self):
        from_dict = self.m.predict(data=self.data, seed=1)
        from_path = self.m.predict(data=str(self.path_test_data), seed=1)
        np.testing.assert_array_equal(from_dict, from_path)

    def test_predictions_are_integer_days_within_window(self):
        pred = self.m.predict(seed=1)
        self.assertTrue(np.issubdtype(pred.dtype, np.integer))
        precursors = np.asarray(self.data["precursor_events"])
        self.assertTrue(np.all(pred >= precursors))
        self.assertTrue(np.all(pred <= max(self.data["N_days"])))

    def test_predict_is_reproducible_with_seed(self):
        np.testing.assert_array_equal(
            self.m.predict(seed=42), self.m.predict(seed=42)
        )

    def test_predict_new_data_without_events(self):
        new = self._first_site_year()
        pred = self.m.predict(data=new, seed=1)
        self.assertEqual(pred.shape, (self.n_draws, new["N"]))
        self.assertTrue(np.all(pred >= np.asarray(new["precursor_events"])))

    def test_predict_does_not_mutate_input(self):
        new = self._first_site_year()
        before = json.dumps(new, sort_keys=True)
        self.m.predict(data=new, seed=1)
        self.assertEqual(json.dumps(new, sort_keys=True), before)

    def test_predictions_track_observed_events(self):
        # Fitted to simulated data, the posterior mean prediction should be
        # close to the observed event days on average.
        pred = self.m.predict(seed=1)
        mae = np.mean(np.abs(pred.mean(axis=0) - np.asarray(self.data["events"])))
        self.assertLess(mae, 15)


if __name__ == "__main__":
    unittest.main()
