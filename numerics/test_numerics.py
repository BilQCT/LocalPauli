"""Fast regressions for errors that previously changed manuscript numerics."""
import itertools
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np
import reproduce as r
import small_core as sc
from export_manuscript import canonical


class NumericalConventions(unittest.TestCase):
    def test_seed_distinguishes_both_modulus_errors(self):
        inputs = np.array(list(itertools.product(range(3), repeat=3)))
        values = sc.gamma2_seed(3)
        for u, expected in [([1, 1, 1], -1), ([2, 2, 0], -1)]:
            index = np.flatnonzero(np.all(inputs == u, axis=1))[0]
            self.assertEqual(values[index], expected)

    def test_graph_matching_respects_decorations(self):
        key = canonical(3, [(0, 1), (1, 2)], [0])
        self.assertEqual(key, canonical(3, [(2, 0), (0, 1)], [2]))
        self.assertNotEqual(key, canonical(3, [(0, 1), (1, 2)], [1]))

    def test_missing_pools_rebuilt_and_corruption_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'data').mkdir()
            with patch.object(r, 'ROOT', root):
                _, _, _, pools = r.pool_small(2, True)
                self.assertEqual(pools['M'].shape, (16, 96))
                self.assertEqual(pools['C'].shape, (16, 432))
                r.pool_small(2, False)
                path = root / 'data/deterministic_2.h5'
                with h5py.File(path, 'r+') as f:
                    f['data'][0, 0] = -1
                with self.assertRaises(AssertionError):
                    r.pool_small(2, False)

    def test_physical_target_normalization_and_gate_order(self):
        digits, _, _ = sc.pauli_data(3)
        for red in ([], [0], [1], [0, 2]):
            actual = r.physical_target(3, [(0, 1), (1, 2)], red)
            expected = sc.resource_coordinates(3, [(0, 1), (1, 2)], red, digits)
            np.testing.assert_allclose(actual, expected, atol=1e-12)
            self.assertAlmostEqual(actual[0], 1)
            self.assertAlmostEqual(actual @ actual, 8)


if __name__ == '__main__':
    unittest.main()
