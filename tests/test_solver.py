from dataclasses import replace

import numpy as np

from mhd_solver.common.eos import primitive_to_conserved
from mhd_solver.hydro1d.exact_sod import exact_sod
from mhd_solver.hydro1d.initial_conditions import sod_shock_tube
from mhd_solver.hydro1d.solver import HydroConfig, apply_outflow_boundaries, solve


def test_outflow_boundaries_copy_edge_cells() -> None:
    conserved = np.arange(15.0).reshape(3, 5)
    extended = apply_outflow_boundaries(conserved)
    np.testing.assert_array_equal(extended[:, :2], conserved[:, :1])
    np.testing.assert_array_equal(extended[:, -2:], conserved[:, -1:])


def test_uniform_state_is_preserved_for_both_orders() -> None:
    primitive = np.array([1.2, 0.3, 0.8])

    def uniform(x: np.ndarray) -> np.ndarray:
        return np.repeat(primitive[:, None], x.size, axis=1)

    for order in (1, 2):
        result = solve(HydroConfig(cells=32, final_time=0.05, order=order), uniform)
        expected = primitive_to_conserved(uniform(result.x), 1.4)
        np.testing.assert_allclose(result.conserved, expected, rtol=1.0e-13, atol=1.0e-13)


def test_sod_solution_remains_physical_and_second_order_is_more_accurate() -> None:
    config = HydroConfig(cells=160, final_time=0.2, cfl=0.7, order=1)
    first = solve(config, sod_shock_tube)
    second = solve(replace(config, order=2), sod_shock_tube)
    exact = exact_sod(second.x, second.time, config.gamma)
    assert first.primitive[0].min() > 0.0 and second.primitive[2].min() > 0.0
    first_density_error = np.mean(np.abs(first.primitive[0] - exact[0]))
    second_density_error = np.mean(np.abs(second.primitive[0] - exact[0]))
    assert second_density_error < first_density_error


def test_closed_uniform_problem_conserves_integrals() -> None:
    def stationary(x: np.ndarray) -> np.ndarray:
        return np.stack((np.ones_like(x), np.zeros_like(x), np.ones_like(x)))

    result = solve(HydroConfig(cells=40, final_time=0.1, order=2), stationary)
    np.testing.assert_allclose(result.final_integrals, result.initial_integrals, atol=1.0e-13)
