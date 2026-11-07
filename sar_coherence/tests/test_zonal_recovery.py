import unittest
import numpy as np
from test_correction import observations, settings
from sarflow.zonal_recovery import classify_series, zonal_classification


def series(pre, event, post):
    return np.array(list(pre) + [event] + list(post), dtype=float)


class ZonalRules(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.obs = observations()

    def classify(self, *columns):
        return classify_series(np.column_stack(columns), self.obs, settings(SOURCE_VALIDITY_CONFIRMED=False))

    def test_decision_tree_classes(self):
        r = self.classify(
            series([.30, .31, .30], .18, [.24, .29, .30]),   # event loss, back to baseline -> 1
            series([.30, .31, .30], .15, [.18, .21, .22]),   # deficit persists, improvement >= .03 -> 2
            series([.30, .31, .30], .15, [.16, .15, .15]),   # no improvement -> 3
            series([.30, .31, .30], .30, [.30, .31, .30]),   # stable -> 4
            series([.30, .31, .30], .20, [.33, .36, .37]),   # event dip then post gain -> 1 (was conflict)
            series([.05, .06, .05], .04, [.05, .05, .06]),   # too low to observe -> 5
        )
        self.assertEqual(r.classes.tolist(), [1, 2, 3, 4, 1, 5])
        self.assertEqual(r.reasons[-1], 4)
        self.assertTrue(r.quality[4] & 8)  # late coherence above baseline + tolerance

    def test_rri_classes_follow_itsa_classes(self):
        r = self.classify(
            series([.30, .31, .30], .18, [.24, .29, .30]),   # ITSA 1, within tolerance -> RRI 4
            series([.30, .31, .30], .20, [.33, .36, .37]),   # ITSA 1, above baseline + tol -> RRI 5
            series([.30, .31, .30], .15, [.18, .21, .22]),   # ITSA 2 -> RRI 3
            series([.30, .31, .30], .15, [.16, .15, .15]),   # ITSA 3, no improvement -> RRI 2
            series([.30, .31, .30], .15, [.14, .12, .12]),   # ITSA 3, below trough -> RRI 1
            series([.30, .31, .30], .30, [.30, .31, .30]),   # ITSA 4 -> RRI 0
            series([.05, .06, .05], .04, [.05, .05, .06]),   # ITSA 5 -> RRI 0
        )
        self.assertEqual(r.classes.tolist(), [1, 1, 2, 3, 3, 4, 5])
        self.assertEqual(r.metrics["rri_class"].tolist(), [4, 5, 3, 2, 1, 0, 0])

    def test_missing_event_is_indeterminate(self):
        y = series([.30, .31, .30], np.nan, [.24, .29, .30])
        r = self.classify(y)
        self.assertEqual(r.classes.tolist(), [5])
        self.assertEqual(r.reasons.tolist(), [2])

    def test_multiscale_falls_back_to_parent_zone(self):
        rng = np.random.default_rng(1)
        T, H, W = len(self.obs), 20, 20
        base = series([.30, .31, .30], .18, [.24, .29, .30])
        stack = np.repeat(base[:, None, None], H, 1).repeat(W, 2).astype('float32')
        stack += rng.normal(0, .005, stack.shape).astype('float32')
        # A 5x5 base zone with only 4 valid pixels must defer to its 10x10 parent.
        stack[:, :5, :5] = np.nan
        stack[:, :2, :2] = base[:, None, None]
        inside = np.ones((H, W), bool)
        owner = np.ones((H, W), 'int64')
        p = zonal_classification(stack, inside, owner, self.obs, settings(), 14.0)
        valid = np.isfinite(stack).any(0)
        self.assertTrue((p.classes[valid] == 1).all())
        self.assertTrue((p.classes[~valid] == 0).all())  # no observations -> NoData
        self.assertTrue((p.zone_size_m[:2, :2] == 140).all())
        self.assertTrue((p.zone_size_m[10:, 10:] == 70).all())
        self.assertFalse((p.classes == 5).any())


if __name__ == '__main__':
    unittest.main()
