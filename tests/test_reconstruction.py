import numpy as np

from mhd_solver.common.reconstruction import minmod


def test_minmod_selects_smallest_same_sign_slope() -> None:
    left = np.array([-2.0, 2.0, -1.0, 1.0])
    right = np.array([-1.0, 3.0, 2.0, -2.0])
    np.testing.assert_allclose(minmod(left, right), [-1.0, 2.0, 0.0, 0.0])

