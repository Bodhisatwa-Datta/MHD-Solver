import unittest

import numpy as np

from mhd_solver.mhd2d.boundaries import fill_boundaries
from mhd_solver.mhd2d.initial_conditions import (
    divergence_perturbation,
    magnetic_rotor,
    orszag_tang_vortex,
)
from mhd_solver.mhd2d.solver import MHD2DConfig, magnetic_divergence, solve
from mhd_solver.reconnection.initial_conditions import harris_current_sheet


class MHD2DSolverTests(unittest.TestCase):
    def test_mixed_boundaries_wrap_x_and_copy_y(self) -> None:
        values = np.arange(12).reshape(1, 3, 4)
        padded = fill_boundaries(values, "periodic", "outflow", ghosts=1)
        np.testing.assert_array_equal(padded[:, 1:-1, 0], values[:, :, -1])
        np.testing.assert_array_equal(padded[:, 1:-1, -1], values[:, :, 0])
        np.testing.assert_array_equal(padded[:, 0, 1:-1], values[:, 0, :])
        np.testing.assert_array_equal(padded[:, -1, 1:-1], values[:, -1, :])

    def test_harris_sheet_is_pressure_balanced_and_divergence_free(self) -> None:
        config = MHD2DConfig(
            nx=24,
            ny=32,
            y_min=-0.5,
            y_max=0.5,
            final_time=0.0,
            boundary_x="periodic",
            boundary_y="outflow",
        )
        dx = (config.x_max - config.x_min) / config.nx
        dy = (config.y_max - config.y_min) / config.ny
        x = config.x_min + (np.arange(config.nx) + 0.5) * dx
        y = config.y_min + (np.arange(config.ny) + 0.5) * dy
        primitive = harris_current_sheet(x, y)
        total_pressure = primitive[4] + 0.5 * primitive[5] ** 2
        np.testing.assert_allclose(total_pressure, 0.7, rtol=0.0, atol=2.0e-15)
        np.testing.assert_allclose(
            primitive[4] / primitive[0], 0.2, rtol=0.0, atol=2.0e-15
        )
        divergence = magnetic_divergence(
            primitive,
            dx,
            dy,
            config.resolved_boundary_x,
            config.resolved_boundary_y,
        )
        self.assertLess(float(np.sqrt(np.mean(divergence**2))), 1.0e-13)

    def test_harris_sheet_short_evolution_remains_physical(self) -> None:
        result = solve(
            MHD2DConfig(
                nx=16,
                ny=32,
                y_min=-0.5,
                y_max=0.5,
                final_time=0.01,
                gamma=5.0 / 3.0,
                cfl=0.3,
                boundary_x="periodic",
                boundary_y="outflow",
                cleaning_speed=2.0,
                cleaning_damping_rate=2.0,
            ),
            harris_current_sheet,
        )
        speed = np.sqrt(result.primitive[1] ** 2 + result.primitive[2] ** 2)
        self.assertGreater(result.primitive[0].min(), 0.0)
        self.assertGreater(result.primitive[4].min(), 0.0)
        self.assertLess(float(speed.max()), 0.02)
        self.assertLess(result.divergence_l2[-1], 1.0e-12)

    def test_magnetic_rotor_initial_state_and_divergence(self) -> None:
        config = MHD2DConfig(nx=32, ny=32, final_time=0.0, gamma=1.4, boundary="outflow")
        dx, dy = 1.0 / config.nx, 1.0 / config.ny
        x, y = (np.arange(config.nx) + 0.5) * dx, (np.arange(config.ny) + 0.5) * dy
        primitive = magnetic_rotor(x, y)
        self.assertGreater(primitive[0].max(), 9.0)
        self.assertAlmostEqual(float(primitive[0].min()), 1.0)
        divergence = magnetic_divergence(primitive, dx, dy, "outflow")
        self.assertLess(float(np.sqrt(np.mean(divergence**2))), 1.0e-13)

    def test_magnetic_rotor_short_evolution_remains_physical(self) -> None:
        result = solve(
            MHD2DConfig(
                nx=16,
                ny=16,
                final_time=0.03,
                gamma=1.4,
                cfl=0.25,
                boundary="outflow",
                cleaning_speed=4.0,
                cleaning_damping_rate=4.0,
            ),
            magnetic_rotor,
        )
        self.assertGreater(result.primitive[0].min(), 0.0)
        self.assertGreater(result.primitive[4].min(), 0.0)
        self.assertTrue(np.all(np.isfinite(result.conserved)))

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
