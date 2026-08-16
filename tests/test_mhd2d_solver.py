import unittest

import numpy as np

from mhd_solver.mhd2d.initial_conditions import divergence_perturbation
from mhd_solver.mhd2d.solver import MHD2DConfig, solve


class MHD2DSolverTests(unittest.TestCase):
    @staticmethod
    def uniform_state(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        shape = (y.size, x.size)
        state = np.array([1.0, 0.2, -0.1, 0.05, 0.8, 0.4, -0.3, 0.2, 0.0])
        return np.broadcast_to(state[:, None, None], (9, *shape)).copy()

    def test_uniform_state_is_preserved_for_both_orders(self) -> None:
        for order in (1, 2):
            config = MHD2DConfig(
                nx=12,
                ny=10,
                final_time=0.02,
                order=order,
                cleaning_damping_rate=0.0,
            )
            result = solve(config, self.uniform_state)
            np.testing.assert_allclose(
                result.primitive,
                self.uniform_state(result.x, result.y),
                rtol=1.0e-13,
                atol=1.0e-13,
            )
            self.assertLess(result.divergence_l2[-1], 1.0e-13)

    def test_glm_damping_reduces_controlled_divergence_error(self) -> None:
        base = dict(nx=16, ny=16, final_time=0.5, cleaning_speed=2.0)
        undamped = solve(
            MHD2DConfig(**base, cleaning_damping_rate=0.0),
            divergence_perturbation,
        )
        damped = solve(
            MHD2DConfig(**base, cleaning_damping_rate=2.0),
            divergence_perturbation,
        )
        self.assertLess(damped.divergence_l2[-1], 0.55 * damped.divergence_l2[0])
        self.assertLess(damped.divergence_l2[-1], undamped.divergence_l2[-1])
        self.assertGreater(damped.primitive[0].min(), 0.0)
        self.assertGreater(damped.primitive[4].min(), 0.0)
