import unittest

import numpy as np

from mhd_solver.mhd2d.equations import (
    conserved_to_primitive,
    fast_magnetosonic_speed,
    hll_flux,
    physical_flux,
    primitive_to_conserved,
)


class MHD2DEquationTests(unittest.TestCase):
    def test_state_round_trip(self) -> None:
        primitive = np.array([1.0, 0.2, -0.1, 0.05, 0.8, 0.4, -0.3, 0.2, 0.1])
        conserved = primitive_to_conserved(primitive, gamma=1.4, cleaning_speed=2.0)
        recovered = conserved_to_primitive(conserved, gamma=1.4, cleaning_speed=2.0)
        np.testing.assert_allclose(recovered, primitive, rtol=1.0e-14, atol=1.0e-14)

    def test_directional_hll_consistency(self) -> None:
        primitive = np.array([1.0, 0.2, -0.1, 0.05, 0.8, 0.4, -0.3, 0.2, 0.1])
        for direction in (0, 1):
            expected = physical_flux(primitive, 1.4, 2.0, direction)
            actual = hll_flux(primitive, primitive, 1.4, 2.0, direction)
            np.testing.assert_allclose(actual, expected, rtol=1.0e-14, atol=1.0e-14)

    def test_zero_field_fast_speed_is_sound_speed(self) -> None:
        primitive = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])
        for direction in (0, 1):
            self.assertAlmostEqual(
                float(fast_magnetosonic_speed(primitive, 1.4, direction)),
                np.sqrt(1.4),
                places=14,
            )
