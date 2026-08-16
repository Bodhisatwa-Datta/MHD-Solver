"""Diagnostics for physical and numerical magnetic diffusion."""

from __future__ import annotations

import numpy as np


def sinusoidal_mode_amplitude(
    field: np.ndarray, x: np.ndarray, wavenumber: float = 2.0 * np.pi
) -> float:
    """Project a cell-centred field onto ``sin(wavenumber*x)``.

    The normalization returns the physical amplitude for a mode sampled on a
    uniform periodic mesh. Any leading dimensions are averaged, allowing both
    one- and two-dimensional field arrays.
    """
    if field.shape[-1] != x.size:
        raise ValueError("the field's final dimension must match x")
    weights = np.sin(wavenumber * x)
    return float(2.0 * np.mean(field * weights))


def effective_resistivity(
    initial_amplitude: float,
    final_amplitude: float,
    elapsed_time: float,
    wavenumber: float = 2.0 * np.pi,
) -> float:
    r"""Infer ``eta`` from modal decay ``A=A_0 exp(-eta*k^2*t)``."""
    if initial_amplitude <= 0.0 or final_amplitude <= 0.0:
        raise ValueError("mode amplitudes must be positive")
    if elapsed_time <= 0.0 or wavenumber == 0.0:
        raise ValueError("elapsed time and wavenumber must be non-zero")
    return float(
        -np.log(final_amplitude / initial_amplitude)
        / (wavenumber**2 * elapsed_time)
    )
