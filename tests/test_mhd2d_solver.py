import unittest

import numpy as np

from mhd_solver.mhd2d.initial_conditions import divergence_perturbation, orszag_tang_vortex
from mhd_solver.mhd2d.solver import MHD2DConfig, magnetic_divergence, solve


class MHD2DSolverTests(unittest.TestCase):
    def test_orszag_tang_initial_field_is_discretely_divergence_free(self) -> None:
        config = MHD2DConfig(nx=32, ny=32, final_time=0.0)
        dx = (config.x_max - config.x_min) / config.nx
        dy = (config.y_max - config.y_min) / config.ny
        x = config.x_min + (np.arange(config.nx) + 0.5) * dx
        y = config.y_min + (np.arange(config.ny) + 0.5) * dy
        primitive = orszag_tang_vortex(x, y)
        divergence = magnetic_divergence(primitive, dx, dy, "periodic")
        self.assertLess(float(np.sqrt(np.mean(divergence**2))), 1.0e-13)

    def test_orszag_tang_short_evolution_remains_physical(self) -> None:
        result = solve(
            MHD2DConfig(nx=16, ny=16, final_time=0.05, cfl=0.3),
            orszag_tang_vortex,
        )
        self.assertGreater(result.primitive[0].min(), 0.0)
        self.assertGreater(result.primitive[4].min(), 0.0)
        self.assertTrue(np.all(np.isfinite(result.conserved)))

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
