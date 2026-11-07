"""Grouped object aggregation must equal per-gid spatial_objects.aggregate exactly.

Small synthetic fixtures only; no research raster is opened. ``segment`` is mocked.
"""
from types import SimpleNamespace
import unittest
from unittest import mock
import numpy as np
import config
from sarflow import study_objects
from sarflow.spatial_objects import aggregate


def settings(**changes):
    c = SimpleNamespace(**{k: v for k, v in vars(config).items() if not k.startswith('__')})
    for k, v in changes.items():
        setattr(c, k, v)
    return c


def reference(stack, ids, c):
    """Original implementation: one full-tile mask scan per object id."""
    keys, counts = np.unique(ids[ids > 0], return_counts=True)
    if not keys.size:
        return keys, counts, np.empty((stack.shape[0], 0))
    return keys, counts, np.column_stack([aggregate(stack, ids == gid, c)[0] for gid in keys])


def assert_same(test, got, want):
    test.assertEqual(got.shape, want.shape)
    test.assertEqual(got.dtype, want.dtype)
    np.testing.assert_array_equal(np.isnan(got), np.isnan(want))
    np.testing.assert_array_equal(got, want)  # NaN==NaN in assert_array_equal; exact values


class GroupedAggregationTests(unittest.TestCase):
    def setUp(self):
        self.c = settings(OBJECT_MIN_VALID_PIXELS=3, OBJECT_MIN_VALID_FRACTION=0.8)

    def check(self, stack, ids, c=None):
        c = c or self.c
        keys, counts, want = reference(stack, ids, c)
        got = study_objects.object_centers(stack, ids, keys, counts, c)
        assert_same(self, got, want)
        return keys, counts, got

    def test_fragmented_singleton_and_varying_sizes(self):
        rng = np.random.default_rng(0)
        ids = np.array([[7, 7, 0, 3, 3, 3],
                        [0, 7, 9, 0, 3, 12],
                        [5, 0, 9, 9, 0, 12],
                        [5, 5, 0, 7, 7, 0]])          # 7 fragmented; 12 two-part; non-contiguous ids
        ids[0, 2] = 42                                 # singleton
        stack = rng.normal(size=(6, *ids.shape)).astype('float32')
        keys, counts, _ = self.check(stack, ids)
        self.assertEqual(keys.tolist(), [3, 5, 7, 9, 12, 42])
        self.assertEqual(counts.tolist(), [4, 3, 5, 3, 2, 1])

    def test_nans_zeros_and_coverage_boundary(self):
        ids = np.array([[1, 1, 1, 1, 1],
                        [2, 2, 2, 2, 2],
                        [3, 3, 0, 4, 0]])
        stack = np.zeros((5, *ids.shape), 'float32')       # valid zeros everywhere
        stack[1, 0, 0] = np.nan                              # 4/5 = 0.8 -> kept, median of zeros
        stack[2, 0, :2] = np.nan                             # 3/5 = 0.6 -> NaN
        stack[3, 1, :] = [0, 0, 5, -1, np.nan]               # median incl. zeros = 0
        stack[4, 1, :] = np.nan                              # all-NaN -> NaN, no warning crash
        stack[1, 2, 0] = np.nan                              # size-2: 1/2 below need -> NaN
        stack[2, 2, 3] = np.nan                              # singleton all-NaN -> NaN
        stack[3, 2, 3] = 0.0                                 # singleton zero is valid
        _, _, got = self.check(stack, ids)
        self.assertEqual(got[1, 0], 0.0)
        self.assertTrue(np.isnan(got[2, 0]))
        self.assertEqual(got[3, 1], 0.0)
        self.assertTrue(np.isnan(got[4, 1]))
        self.assertTrue(np.isnan(got[1, 2]))
        self.assertTrue(np.isnan(got[2, 3]))
        self.assertEqual(got[3, 3], 0.0)

    def test_min_valid_pixels_capped_at_group_size(self):
        ids = np.array([[1, 2, 2, 0, 3, 3, 3, 3]])
        stack = np.arange(3 * ids.size, dtype='float64').reshape(3, *ids.shape)
        stack[1, 0, 4] = np.nan                       # size 4: 3 valid, fraction 0.75 passes 0.5
        c = settings(OBJECT_MIN_VALID_PIXELS=4, OBJECT_MIN_VALID_FRACTION=0.5)
        _, _, got = self.check(stack, ids, c)         # need = min(4, size) for sizes 1, 2
        self.assertTrue(np.isfinite(got[:, :2]).all())
        self.assertTrue(np.isnan(got[1, 2]))          # 3 valid < need 4

    def test_random_property_with_chunking(self):
        rng = np.random.default_rng(1)
        for trial in range(25):
            shape = (int(rng.integers(1, 12)), int(rng.integers(1, 12)))
            ids = rng.integers(0, int(rng.integers(1, 9)), size=shape) * int(rng.integers(1, 4))
            stack = rng.normal(size=(int(rng.integers(1, 7)), *shape)).astype('float32')
            stack[rng.random(stack.shape) < 0.25] = np.nan
            stack[rng.random(stack.shape) < 0.15] = 0.0
            budget = int(rng.choice([1, 7, 1 << 24]))   # forces per-object and multi-object chunks
            with mock.patch.object(study_objects, '_GATHER_BUDGET', budget):
                with self.subTest(trial=trial, budget=budget):
                    self.check(stack, ids)

    def test_non_contiguous_stack(self):
        rng = np.random.default_rng(2)
        ids = rng.integers(0, 4, size=(6, 5))
        stack = rng.normal(size=(4, 5, 6))[:, ::-1, :].transpose(0, 2, 1)   # (4, 6, 5) view
        self.assertFalse(stack.flags.c_contiguous)
        self.check(stack, ids)


class ObjectSeriesTests(unittest.TestCase):
    def run_series(self, stack, ids, c):
        obs = [SimpleNamespace(phase=p) for p in ('pre', 'pre', 'post', 'post')][:stack.shape[0]]
        inside = np.ones(ids.shape, bool)
        with mock.patch.object(study_objects, 'cell_sizes', return_value=(1, None)), \
             mock.patch.object(study_objects, 'pre_change_only', return_value=None), \
             mock.patch.object(study_objects, 'day_offsets', return_value=np.zeros(len(obs))), \
             mock.patch.object(study_objects, 'segment', return_value=ids):
            return study_objects.object_series(stack, inside, obs, c, 1.0, 0, 0, 1)

    def test_matches_original_outputs(self):
        c = settings(OBJECT_MIN_VALID_PIXELS=3, OBJECT_MIN_VALID_FRACTION=0.8)
        rng = np.random.default_rng(3)
        ids = np.array([[4, 4, 0, 9], [4, 0, 9, 9], [0, 2, 2, 11]])
        stack = rng.normal(size=(4, *ids.shape)).astype('float32')
        stack[2, 0, 0] = np.nan
        stack[3, 1, 2] = 0.0
        out_ids, keys, counts, series = self.run_series(stack, ids, c)
        want_keys, want_counts, want = reference(stack, ids, c)
        self.assertIs(out_ids, ids)
        np.testing.assert_array_equal(keys, want_keys)
        np.testing.assert_array_equal(counts, want_counts)
        assert_same(self, series, want)

    def test_no_objects(self):
        c = settings()
        ids = np.zeros((3, 3), int)
        stack = np.ones((4, 3, 3), 'float32')
        _, keys, counts, series = self.run_series(stack, ids, c)
        self.assertEqual(keys.size, 0)
        self.assertEqual(counts.size, 0)
        self.assertEqual(series.shape, (4, 0))


if __name__ == '__main__':
    unittest.main()
