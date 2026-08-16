import numpy as np
import unittest

from mhd_solver.common.reconstruction import minmod


class ReconstructionTests(unittest.TestCase):
    def test_minmod_selects_smallest_same_sign_slope(self) -> None:
        left = np.array([-2.0, 2.0, -1.0, 1.0])
        right = np.array([-1.0, 3.0, 2.0, -2.0])
        np.testing.assert_allclose(minmod(left, right), [-1.0, 2.0, 0.0, 0.0])
