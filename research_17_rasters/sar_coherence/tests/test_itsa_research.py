"""Controlled software tests only; never substitutes for real input observations.

Written for a separately authorized validation run. Not executed during revision.
"""
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import config
from sarflow.temporal import Observation, configured_pair, day_offsets
from sarflow.itsa import fit_itsa
from sarflow.recovery import assess
from sarflow.spatial_objects import (aggregate, segment, pre_change_only,
                                     region_assessments, regularize, cell_sizes)


def settings(**changes):
    result = SimpleNamespace(**{k: v for k, v in vars(config).items() if not k.startswith("__")})
    result.EVENT_ORIGIN = "2022-09-25"
    result.SOURCE_VALIDITY_CONFIRMED = True
    result.SOURCE_VALIDITY_NOTES = "Synthetic unit-test magnitudes only"
    for name, value in changes.items():
        setattr(result, name, value)
    return result


def observations():
    origin = pd.Timestamp("2022-09-25")
    phases = ["pre"]*8 + ["event"] + ["post"]*8
    times = [-90, -78, -66, -54, -42, -30, -18, -6, 6, 18, 30, 42, 54, 66, 78, 90, 102]
    result = []
    for i, (phase, day) in enumerate(zip(phases, times)):
        end = origin + pd.Timedelta(days=day)
        start = end-pd.Timedelta(days=12)
        result.append(Observation(f"o{i}", Path(f"unit_{i}.tif"), phase, start, end, start, end, end))
    return result


def sequence(slope=.003, loss=.20, event=.4):
    obs = observations()
    t = day_offsets(obs, "2022-09-25")
    y = np.full(17, .6)
    y[8] = event
    y[9:] = .6-loss+slope*t[9:]
    return y


class ITSARecoveryTests(unittest.TestCase):
    def setUp(self):
        self.c, self.obs = settings(), observations()

    def result(self, y):
        return assess(y[:, None], self.obs, self.c)

    def test_complete_recovery(self):
        self.assertEqual(int(self.result(sequence()).classes[0]), 1)

    def test_partial_recovery(self):
        self.assertEqual(int(self.result(sequence(slope=.0008)).classes[0]), 2)

    def test_no_recovery(self):
        self.assertEqual(int(self.result(sequence(slope=0)).classes[0]), 3)

    def test_further_degradation(self):
        self.assertEqual(int(self.result(sequence(slope=-.0005)).classes[0]), 3)

    def test_no_loss(self):
        self.assertEqual(int(self.result(sequence(slope=0, loss=0, event=.6)).classes[0]), 4)

    def test_contradictory_event_is_not_no_loss(self):
        result = self.result(sequence(event=.8))
        self.assertEqual(int(result.classes[0]), 5)
        self.assertEqual(int(result.reasons[0]), 5)

    def test_wide_uncertainty_is_not_automatic_rejection(self):
        self.c.UNCERTAINTY_MULTIPLIER = 1e9
        y = sequence()
        y += np.array([.003, -.003]*8+[.003])
        result = self.result(y)
        self.assertEqual(int(result.classes[0]), 1)
        self.assertTrue(int(result.quality[0]) & 1)

    def test_missing_support(self):
        y = sequence()
        y[:4] = np.nan
        self.assertEqual(int(self.result(y).classes[0]), 5)

    def test_missing_latest_observation_breaks_confirmation(self):
        y = sequence()
        y[-1] = np.nan
        self.assertEqual(int(self.result(y).classes[0]), 5)

    def test_middle_missing_is_not_interpolated(self):
        y = sequence()
        y[12] = np.nan
        result = self.result(y)
        self.assertEqual(int(result.classes[0]), 1)
        self.assertEqual(int(result.metrics["post_count"][0]), 7)

    def test_model_dependence(self):
        # Identical PRE/event; changing the POST fitted trajectory changes the outcome.
        self.assertNotEqual(self.result(sequence()).classes[0], self.result(sequence(slope=0)).classes[0])

    def test_against_independent_closed_form_segments(self):
        y = sequence()
        t = day_offsets(self.obs, self.c.EVENT_ORIGIN)
        pre = np.polyfit(t[:8], y[:8], 1)
        post = np.polyfit(t[9:], y[9:], 1)
        expected = [pre[1], pre[0], post[1]-pre[1], post[0]-pre[0]]
        np.testing.assert_allclose(fit_itsa(y[:, None], self.obs, self.c).beta[:, 0], expected, atol=1e-10)

    def test_event_is_supporting_not_regression_point(self):
        a, b = sequence(), sequence(event=.2)
        np.testing.assert_allclose(fit_itsa(a[:, None], self.obs, self.c).beta,
                                   fit_itsa(b[:, None], self.obs, self.c).beta)

    def test_all_eight_pre_dates_are_used(self):
        y = sequence()
        before = fit_itsa(y[:, None], self.obs, self.c).beta.copy()
        y[6] += .06
        after = fit_itsa(y[:, None], self.obs, self.c).beta
        self.assertFalse(np.allclose(before, after))

    def test_zero_is_observed_not_structurally_missing(self):
        result = self.result(np.zeros(17))
        self.assertEqual(int(result.classes[0]), 5)
        self.assertEqual(int(result.metrics["pre_count"][0]), 8)

    def test_structural_nodata(self):
        self.assertEqual(int(self.result(np.full(17, np.nan)).classes[0]), 0)

    def test_supported_crossing_only(self):
        result = self.result(sequence())
        crossing = float(result.metrics["crossing_days"][0])
        self.assertTrue(18 <= crossing <= 102)
        self.assertTrue(np.isnan(self.result(sequence(slope=0)).metrics["crossing_days"][0]))

    def test_filename_pair_parse(self):
        c = SimpleNamespace(PAIR_DATES={"Pre1": ("", "")})
        a, b = configured_pair(c, "Pre1", Path("coherence_20220721_20220802.tif"))
        self.assertEqual((str(a.date()), str(b.date())), ("2022-07-21", "2022-08-02"))


class ObjectTests(unittest.TestCase):
    def setUp(self):
        self.c, self.obs = settings(), observations()
        self.stack = np.broadcast_to(sequence()[:, None, None], (17, 6, 6)).copy()
        self.mask = np.ones((6, 6), bool)
        self.barriers = np.ones((6, 6), dtype="uint8")

    def ids(self, row=0, col=0, barriers=None):
        pre = self.stack[:8]
        change = pre_change_only(pre, np.arange(8))
        return segment(pre, None, change, self.mask, self.barriers if barriers is None else barriers,
                       6, row, col, 100, self.c)

    def test_contiguous_homogeneous_object(self):
        self.assertEqual(len(np.unique(self.ids())), 1)

    def test_pre_heterogeneity_splits_objects(self):
        self.stack[:8, :, 3:] = .8
        self.assertEqual(len(np.unique(self.ids())), 2)

    def test_nodata_hole_never_filled(self):
        self.mask[2, 2] = False
        self.assertEqual(int(self.ids()[2, 2]), 0)

    def test_external_barrier(self):
        self.barriers[:, 3:] = 2
        self.assertEqual(len(np.unique(self.ids())), 2)

    def test_object_valid_fraction(self):
        self.stack[0, :3] = np.nan
        series, count, fraction, _ = aggregate(self.stack, self.mask, self.c)
        self.assertTrue(np.isnan(series[0]))
        self.assertEqual(int(count[0]), 18)
        self.assertEqual(float(fraction[0]), .5)

    def test_region_itsa_recovered(self):
        ids = self.ids()
        records = region_assessments(self.stack, ids, self.obs, self.c, 200.)
        self.assertEqual(int(next(iter(records.values()))[0].classes[0]), 1)

    def test_small_supported_object_is_flagged_not_removed(self):
        ids = self.ids()
        records = region_assessments(self.stack, ids, self.obs, self.c, 1.)
        result = next(iter(records.values()))[0]
        self.assertEqual(int(result.classes[0]), 1)
        self.assertTrue(int(result.quality[0]) & 4)

    def test_global_ids_do_not_collide_across_domains(self):
        self.assertFalse(set(np.unique(self.ids())).intersection(np.unique(self.ids(row=6))))

    def test_domain_geometry_independent_of_io_window_size(self):
        a = cell_sizes(200., self.c)
        self.c.WINDOW_SIZE = 512
        self.assertEqual(a, cell_sizes(200., self.c))

    def test_high_confidence_small_degraded_patch_is_preserved(self):
        ids = self.ids()
        ids[:, 3:] = 4
        self.stack[9:, :, 3:] = .4
        records = region_assessments(self.stack, ids, self.obs, self.c, 200.)
        for entry in records.values():
            entry[0].quality[:] = 0
        after, _ = regularize(self.stack, ids, records, self.barriers, self.obs, self.c, 200.)
        np.testing.assert_array_equal(after, ids)

    def test_indeterminate_never_promoted_by_neighbors(self):
        ids = self.ids()
        records = region_assessments(self.stack, ids, self.obs, self.c, 200.)
        next(iter(records.values()))[0].classes[:] = 5
        after, results = regularize(self.stack, ids, records, self.barriers, self.obs, self.c, 200.)
        self.assertEqual(int(next(iter(results.values()))[0].classes[0]), 5)
        np.testing.assert_array_equal(after, ids)


if __name__ == "__main__":
    unittest.main()
