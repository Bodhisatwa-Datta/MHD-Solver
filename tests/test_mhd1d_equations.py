import unittest

import numpy as np

from mhd_solver.common.eos import UnphysicalStateError
from mhd_solver.mhd1d.equations import (
    conserved_to_primitive,
    fast_magnetosonic_speed,
    hll_flux,
    physical_flux,
    primitive_to_conserved,
)


class MHD1DEquationTests(unittest.TestCase):
    def test_primitive_conserved_round_trip(self) -> None:
        primitive = np.array(
            [[1.0, 0.4], [0.2, -0.1], [0.1, 0.3], [0.0, -0.2],
             [1.0, 0.3], [0.5, -0.2], [0.1, 0.4]]
        )
        conserved = primitive_to_conserved(primitive, gamma=1.4, longitudinal_field=0.75)
        recovered = conserved_to_primitive(conserved, gamma=1.4, longitudinal_field=0.75)
        np.testing.assert_allclose(recovered, primitive, rtol=1.0e-14, atol=1.0e-14)

    def test_zero_field_fast_speed_reduces_to_sound_speed(self) -> None:
        primitive = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0])
        speed = fast_magnetosonic_speed(primitive, gamma=1.4, longitudinal_field=0.0)
        self.assertAlmostEqual(float(speed), np.sqrt(1.4), places=14)

    def test_hll_flux_is_consistent_for_equal_states(self) -> None:
        primitive = np.array([1.0, 0.2, -0.1, 0.05, 1.0, 0.4, -0.2])
        expected = physical_flux(primitive, gamma=1.4, longitudinal_field=0.75)
        actual = hll_flux(primitive, primitive, gamma=1.4, longitudinal_field=0.75)
        np.testing.assert_allclose(actual, expected, rtol=1.0e-14, atol=1.0e-14)

    def test_negative_gas_pressure_is_rejected(self) -> None:
        primitive = np.array([1.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0])
        with self.assertRaisesRegex(UnphysicalStateError, "pressure"):
            primitive_to_conserved(primitive, gamma=1.4, longitudinal_field=0.75)
