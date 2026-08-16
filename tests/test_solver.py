from dataclasses import replace
import unittest

import numpy as np

from mhd_solver.common.eos import primitive_to_conserved
from mhd_solver.hydro1d.exact_sod import exact_sod
from mhd_solver.hydro1d.initial_conditions import smooth_density_wave, sod_shock_tube
from mhd_solver.hydro1d.solver import (
    HydroConfig,
    apply_outflow_boundaries,
    apply_periodic_boundaries,
    solve,
)


class SolverTests(unittest.TestCase):
    def test_exact_sod_star_state_matches_standard_values(self) -> None:
        sample = exact_sod(np.array([0.6, 0.75]), 0.2)
        np.testing.assert_allclose(sample[1], [0.92745262, 0.92745262], rtol=1.0e-7)
        np.testing.assert_allclose(sample[2], [0.30313018, 0.30313018], rtol=1.0e-7)

    def test_outflow_boundaries_copy_edge_cells(self) -> None:
        conserved = np.arange(15.0).reshape(3, 5)
        extended = apply_outflow_boundaries(conserved)
        np.testing.assert_array_equal(extended[:, :2], np.repeat(conserved[:, :1], 2, axis=1))
        np.testing.assert_array_equal(extended[:, -2:], np.repeat(conserved[:, -1:], 2, axis=1))

    def test_periodic_boundaries_wrap_opposite_edges(self) -> None:
        conserved = np.arange(15.0).reshape(3, 5)
        extended = apply_periodic_boundaries(conserved)
        np.testing.assert_array_equal(extended[:, :2], conserved[:, -2:])
        np.testing.assert_array_equal(extended[:, -2:], conserved[:, :2])

    def test_uniform_state_is_preserved_for_both_orders(self) -> None:
        primitive = np.array([1.2, 0.3, 0.8])

        def uniform(x: np.ndarray) -> np.ndarray:
            return np.repeat(primitive[:, None], x.size, axis=1)

        for order in (1, 2):
            result = solve(HydroConfig(cells=32, final_time=0.05, order=order), uniform)
            expected = primitive_to_conserved(uniform(result.x), 1.4)
            np.testing.assert_allclose(result.conserved, expected, rtol=1.0e-13, atol=1.0e-13)

    def test_sod_solution_is_physical_and_second_order_more_accurate(self) -> None:
        config = HydroConfig(cells=160, final_time=0.2, cfl=0.7, order=1)
        first = solve(config, sod_shock_tube)
        second = solve(replace(config, order=2), sod_shock_tube)
        exact = exact_sod(second.x, second.time, config.gamma)
        self.assertGreater(first.primitive[0].min(), 0.0)
        self.assertGreater(second.primitive[2].min(), 0.0)
        first_error = np.mean(np.abs(first.primitive[0] - exact[0]))
        second_error = np.mean(np.abs(second.primitive[0] - exact[0]))
        self.assertLess(second_error, first_error)

    def test_closed_uniform_problem_conserves_integrals(self) -> None:
        def stationary(x: np.ndarray) -> np.ndarray:
            return np.stack((np.ones_like(x), np.zeros_like(x), np.ones_like(x)))

        result = solve(HydroConfig(cells=40, final_time=0.1, order=2), stationary)
        np.testing.assert_allclose(result.final_integrals, result.initial_integrals, atol=1.0e-13)

    def test_periodic_density_wave_converges_under_refinement(self) -> None:
        errors = []
        for cells in (40, 80):
            config = HydroConfig(
                cells=cells, final_time=0.25, cfl=0.6, order=2, boundary="periodic"
            )
            result = solve(config, smooth_density_wave)
            exact = smooth_density_wave(result.x, time=result.time)
            errors.append(np.mean(np.abs(result.primitive[0] - exact[0])))
            np.testing.assert_allclose(
                result.final_integrals, result.initial_integrals, atol=2.0e-13
            )
        self.assertGreater(errors[0] / errors[1], 3.0)
