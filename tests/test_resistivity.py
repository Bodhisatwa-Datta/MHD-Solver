import unittest

import numpy as np

from mhd_solver.mhd2d.solver import MHD2DConfig, solve
from mhd_solver.resistive.initial_conditions import magnetic_diffusion_wave
from mhd_solver.resistive.resistivity import resistive_rhs


class ResistivityTests(unittest.TestCase):
    def test_sinusoidal_diffusion_operator_is_second_order(self) -> None:
        errors = []
        for cells in (32, 64):
            x = (np.arange(cells) + 0.5) / cells
            y = (np.arange(cells) + 0.5) / cells
            primitive = magnetic_diffusion_wave(x, y)
            rhs = resistive_rhs(primitive, 0.1, 1.0 / cells, 1.0 / cells, "periodic")
            exact = -0.1 * (2.0 * np.pi) ** 2 * primitive[6]
            errors.append(float(np.mean(np.abs(rhs[6] - exact))))
            self.assertLess(abs(float(np.sum(rhs[4]))), 1.0e-14)
        self.assertGreater(errors[0] / errors[1], 3.9)

    def test_integrated_diffusion_converges_and_conserves_total_energy(self) -> None:
        errors = []
        for cells in (16, 32):
            config = MHD2DConfig(
                nx=cells,
                ny=cells,
                final_time=0.02,
                gamma=5.0 / 3.0,
                cfl=0.35,
                boundary="periodic",
                cleaning_speed=2.0,
                cleaning_damping_rate=0.0,
                resistivity=0.1,
            )
            result = solve(config, magnetic_diffusion_wave)
            exact = magnetic_diffusion_wave(
                result.x, result.y, time=result.time, resistivity=config.resistivity
            )
            errors.append(float(np.mean(np.abs(result.primitive[6] - exact[6]))))
            self.assertLess(
                abs(float(result.final_integrals[4] - result.initial_integrals[4])),
                2.0e-13,
            )
        self.assertGreater(errors[0] / errors[1], 3.0)
