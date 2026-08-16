import numpy as np
import pytest

from mhd_solver.common.eos import (
    UnphysicalStateError,
    conserved_to_primitive,
    primitive_to_conserved,
)
from mhd_solver.common.riemann import euler_flux, hll_flux


def test_primitive_conserved_round_trip() -> None:
    primitive = np.array([[1.0, 0.5], [0.25, -0.1], [1.0, 0.2]])
    recovered = conserved_to_primitive(primitive_to_conserved(primitive, 1.4), 1.4)
    np.testing.assert_allclose(recovered, primitive, rtol=1.0e-14, atol=1.0e-14)


def test_negative_pressure_is_rejected() -> None:
    with pytest.raises(UnphysicalStateError, match="pressure"):
        primitive_to_conserved(np.array([1.0, 0.0, -1.0]), 1.4)


def test_stationary_uniform_flux() -> None:
    state = np.array([1.0, 0.0, 1.0])
    np.testing.assert_allclose(euler_flux(state, 1.4), [0.0, 1.0, 0.0])
    np.testing.assert_allclose(hll_flux(state, state, 1.4), euler_flux(state, 1.4))

