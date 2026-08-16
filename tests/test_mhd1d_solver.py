import unittest

import numpy as np

from mhd_solver.mhd1d.equations import primitive_to_conserved
from mhd_solver.mhd1d.initial_conditions import (
    brio_wu_shock_tube,
    circularly_polarized_alfven_wave,
)
from mhd_solver.mhd1d.solver import MHD1DConfig, solve


class MHD1DSolverTests(unittest.TestCase):
    def test_alfven_wave_returns_after_one_period(self) -> None:
        x = np.linspace(0.0, 1.0, 17, endpoint=False)
        initial = circularly_polarized_alfven_wave(x, time=0.0)
        one_period = circularly_polarized_alfven_wave(x, time=1.0)
        np.testing.assert_allclose(one_period, initial, atol=1.0e-14)

    def test_uniform_state_is_preserved_for_both_orders(self) -> None:
        state = np.array([1.0, 0.2, -0.1, 0.05, 1.0, 0.4, -0.2])

        def uniform(x: np.ndarray) -> np.ndarray:
            return np.repeat(state[:, None], x.size, axis=1)

        for order in (1, 2):
            config = MHD1DConfig(
                cells=32,
                final_time=0.02,
                gamma=1.4,
                order=order,
                longitudinal_field=0.75,
            )
            result = solve(config, uniform)
            expected = primitive_to_conserved(
                uniform(result.x), config.gamma, config.longitudinal_field
            )
            np.testing.assert_allclose(result.conserved, expected, atol=2.0e-14)

    def test_brio_wu_remains_physical(self) -> None:
        config = MHD1DConfig(cells=200, final_time=0.1, cfl=0.4, order=2)
        result = solve(config, brio_wu_shock_tube)
        self.assertGreater(result.primitive[0].min(), 0.0)
        self.assertGreater(result.primitive[4].min(), 0.0)
        self.assertTrue(np.all(np.isfinite(result.conserved)))

    def test_periodic_uniform_state_conserves_all_evolved_variables(self) -> None:
        state = np.array([0.8, -0.1, 0.2, 0.0, 0.7, -0.3, 0.15])

        def uniform(x: np.ndarray) -> np.ndarray:
            return np.repeat(state[:, None], x.size, axis=1)

        result = solve(
            MHD1DConfig(
                cells=40,
                final_time=0.03,
                gamma=1.4,
                boundary="periodic",
            ),
            uniform,
        )
        np.testing.assert_allclose(
            result.final_integrals, result.initial_integrals, atol=2.0e-14
        )

    def test_alfven_wave_second_order_error_falls_under_refinement(self) -> None:
        errors = []
        for cells in (40, 80):
            config = MHD1DConfig(
                cells=cells,
                x_min=0.0,
                x_max=1.0,
                final_time=0.25,
                gamma=5.0 / 3.0,
                cfl=0.4,
                order=2,
                boundary="periodic",
                longitudinal_field=1.0,
            )
            result = solve(config, circularly_polarized_alfven_wave)
            exact = circularly_polarized_alfven_wave(result.x, time=result.time)
            field_error = np.sqrt(
                (result.primitive[5] - exact[5]) ** 2
                + (result.primitive[6] - exact[6]) ** 2
            )
            errors.append(float(np.mean(field_error)))
            np.testing.assert_allclose(
                result.final_integrals, result.initial_integrals, atol=3.0e-13
            )
        self.assertGreater(errors[0] / errors[1], 3.0)
