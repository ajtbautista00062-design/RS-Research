"""Independent seven-observation fixtures for the explicit observability amendment."""
from pathlib import Path
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

import config
from sarflow.recovery import assess
from sarflow.temporal import Observation, day_offsets


def settings(**changes):
    values = {k: v for k, v in vars(config).items() if k.isupper()}
    values.update(EVENT_ORIGIN="2022-09-25T09:30:00", ITSA_MIN_PRE=3,
                  ITSA_MIN_POST=3, BASELINE_ELIGIBILITY_MODE="signal_evidence",
                  OBSERVABILITY_NOISE_MULTIPLIER=3.0,
                  SOURCE_VALIDITY_CONFIRMED=True)
    values.update(changes)
    return SimpleNamespace(**values)


def observations():
    pairs = [("20220721", "20220802"), ("20220802", "20220814"),
             ("20220907", "20220919"), ("20220919", "20221001"),
             ("20221001", "20221013"), ("20221013", "20221025"),
             ("20221025", "20221106")]
    phases = ["pre"] * 3 + ["event"] + ["post"] * 3
    result = []
    for i, ((left, right), phase) in enumerate(zip(pairs, phases)):
        a, b = pd.Timestamp(left), pd.Timestamp(right)
        result.append(Observation(str(i), Path(f"fixture_{i}.tif"), phase,
                                  a, b, a, b, b))
    return result


def series(baseline=.3, loss=.2, slope=0., event=None):
    days = day_offsets(observations(), "2022-09-25T09:30:00")
    result = np.full(7, baseline, dtype=float)
    result[3] = baseline - loss if event is None else event
    result[4:] = baseline - loss + slope * days[4:]
    return result


class CorrectionTests(unittest.TestCase):
    def evaluate(self, values, **changes):
        return assess(np.asarray(values)[:, None], observations(), settings(**changes))

    def test_low_baseline_persistent_loss_with_signal(self):
        result = self.evaluate(series())
        self.assertEqual(int(result.classes[0]), 3)
        self.assertTrue(int(result.quality[0]) & 64)
        self.assertFalse(result.gates["baseline_absolute_pass"][0])
        self.assertTrue(result.gates["baseline_signal_pass"][0])
        self.assertTrue(result.gates["observability_noise_pass"][0])

    def test_low_baseline_insufficient_dynamic_range(self):
        result = self.evaluate(np.full(7, .02))
        self.assertEqual(int(result.classes[0]), 5)
        self.assertEqual(int(result.reasons[0]), 4)

    def test_low_baseline_noisy_uninformative_series(self):
        result = self.evaluate([.01, .35, .02, .05, .02, .25, .01])
        self.assertEqual(int(result.classes[0]), 5)
        self.assertFalse(result.gates["observability_noise_pass"][0])

    def test_stable_low_baseline_no_loss(self):
        self.assertEqual(int(self.evaluate(np.full(7, .3)).classes[0]), 4)

    def test_partial_recovery(self):
        self.assertEqual(int(self.evaluate(series(slope=.003)).classes[0]), 2)

    def test_sustained_recovery(self):
        self.assertEqual(int(self.evaluate(series(slope=.0052)).classes[0]), 1)

    def test_transient_event_dip_alone_stays_unresolved(self):
        result = self.evaluate(series(loss=-.05, event=.1))
        self.assertEqual(int(result.classes[0]), 5)
        self.assertEqual(int(result.reasons[0]), 5)
        self.assertTrue(result.gates["transient_candidate"][0])

    def test_opposite_event_and_persistent_loss_stays_conflicting(self):
        result = self.evaluate(series(baseline=.6, event=.8))
        self.assertEqual(int(result.classes[0]), 5)
        self.assertEqual(int(result.reasons[0]), 5)

    def test_recovery_ambiguity_has_specific_reason(self):
        result = self.evaluate(series(baseline=.6, slope=.00125))
        self.assertEqual(int(result.classes[0]), 5)
        self.assertEqual(int(result.reasons[0]), 9)
        self.assertTrue(result.gates["recovery_ready"][0])
        self.assertTrue(result.gates["ambiguous_recovery"][0])

    def test_near_disturbance_threshold_remains_reason_six(self):
        result = self.evaluate(series(baseline=.6, loss=.03))
        self.assertEqual(int(result.classes[0]), 5)
        self.assertEqual(int(result.reasons[0]), 6)

    def test_missing_observation_is_not_interpolated(self):
        y = series()
        y[1] = np.nan
        result = self.evaluate(y)
        self.assertEqual(int(result.classes[0]), 5)
        self.assertEqual(int(result.reasons[0]), 2)
        self.assertEqual(int(result.metrics["pre_count"][0]), 2)

    def test_source_nodata_stays_nodata(self):
        result = self.evaluate(np.full(7, np.nan))
        self.assertEqual(int(result.classes[0]), 0)
        self.assertEqual(int(result.reasons[0]), 1)

    def test_legitimate_zero_counts_as_observation(self):
        result = self.evaluate(np.zeros(7))
        self.assertEqual(int(result.classes[0]), 5)
        self.assertEqual(int(result.metrics["pre_count"][0]), 3)
        self.assertEqual(int(result.metrics["post_count"][0]), 3)

    def test_historical_configuration_stays_absolute(self):
        result = self.evaluate(series(), BASELINE_ELIGIBILITY_MODE="legacy_absolute")
        self.assertEqual(int(result.classes[0]), 5)
        self.assertEqual(int(result.reasons[0]), 4)

    def test_missing_mode_setting_retains_historical_screen(self):
        c = settings()
        del c.BASELINE_ELIGIBILITY_MODE
        result = assess(series()[:, None], observations(), c)
        self.assertEqual(int(result.classes[0]), 5)
        self.assertEqual(int(result.reasons[0]), 4)

    def test_event_stays_out_of_ols_and_actual_days_are_used(self):
        y = series(slope=.003)
        a = self.evaluate(y)
        y[3] = .01
        b = self.evaluate(y)
        np.testing.assert_allclose(a.fit.beta, b.fit.beta, atol=1e-12)
        np.testing.assert_allclose(a.fit.beta[:, 0], [.3, 0., -.2, .003], atol=1e-12)
        self.assertEqual(int(a.fit.df[0]), 2)


if __name__ == "__main__":
    unittest.main()
