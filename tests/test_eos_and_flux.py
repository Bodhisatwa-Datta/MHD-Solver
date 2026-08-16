import unittest

import numpy as np

from mhd_solver.common.eos import (
    UnphysicalStateError,
    conserved_to_primitive,
    primitive_to_conserved,
)
from mhd_solver.common.riemann import euler_flux, hll_flux


class EquationOfStateAndFluxTests(unittest.TestCase):
    def test_primitive_conserved_round_trip(self) -> None:
        primitive = np.array([[1.0, 0.5], [0.25, -0.1], [1.0, 0.2]])
        recovered = conserved_to_primitive(primitive_to_conserved(primitive, 1.4), 1.4)
        np.testing.assert_allclose(recovered, primitive, rtol=1.0e-14, atol=1.0e-14)

    def test_negative_pressure_is_rejected(self) -> None:
        with self.assertRaisesRegex(UnphysicalStateError, "pressure"):
            primitive_to_conserved(np.array([1.0, 0.0, -1.0]), 1.4)

    def test_stationary_uniform_flux(self) -> None:
        state = np.array([1.0, 0.0, 1.0])
        np.testing.assert_allclose(euler_flux(state, 1.4), [0.0, 1.0, 0.0])
        np.testing.assert_allclose(hll_flux(state, state, 1.4), euler_flux(state, 1.4))
