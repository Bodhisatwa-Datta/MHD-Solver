"""Exact self-similar Riemann solution for the standard Sod problem."""

from __future__ import annotations

import numpy as np


def _pressure_function(
    pressure: float, density: float, initial_pressure: float, sound: float, gamma: float
) -> tuple[float, float]:
    if pressure > initial_pressure:
        coefficient_a = 2.0 / ((gamma + 1.0) * density)
        coefficient_b = (gamma - 1.0) / (gamma + 1.0) * initial_pressure
        root = np.sqrt(coefficient_a / (pressure + coefficient_b))
        value = (pressure - initial_pressure) * root
        derivative = root * (1.0 - 0.5 * (pressure - initial_pressure) / (pressure + coefficient_b))
        return value, derivative
    exponent = (gamma - 1.0) / (2.0 * gamma)
    ratio = pressure / initial_pressure
    value = 2.0 * sound / (gamma - 1.0) * (ratio**exponent - 1.0)
    derivative = (1.0 / (density * sound)) * ratio ** (
        -(gamma + 1.0) / (2.0 * gamma)
    )
    return value, derivative


def _star_state(gamma: float) -> tuple[float, float]:
    rho_l, p_l, rho_r, p_r = 1.0, 1.0, 0.125, 0.1
    sound_l = np.sqrt(gamma * p_l / rho_l)
    sound_r = np.sqrt(gamma * p_r / rho_r)
    pressure = 0.5 * (p_l + p_r)
    for _ in range(50):
        left, dleft = _pressure_function(pressure, rho_l, p_l, sound_l, gamma)
        right, dright = _pressure_function(pressure, rho_r, p_r, sound_r, gamma)
        updated = max(pressure - (left + right) / (dleft + dright), 1.0e-12)
        if abs(updated - pressure) <= 1.0e-12 * (updated + pressure):
            pressure = updated
            break
        pressure = updated
    left, _ = _pressure_function(pressure, rho_l, p_l, sound_l, gamma)
    right, _ = _pressure_function(pressure, rho_r, p_r, sound_r, gamma)
    return pressure, 0.5 * (right - left)


def exact_sod(
    x: np.ndarray, time: float, gamma: float = 1.4, discontinuity: float = 0.5
) -> np.ndarray:
    """Sample the exact standard Sod solution as ``(rho, velocity, pressure)``."""
    if time <= 0.0:
        from mhd_solver.hydro1d.initial_conditions import sod_shock_tube

        return sod_shock_tube(x, discontinuity)

    rho_l, p_l, rho_r, p_r = 1.0, 1.0, 0.125, 0.1
    sound_l = np.sqrt(gamma * p_l / rho_l)
    sound_r = np.sqrt(gamma * p_r / rho_r)
    p_star, velocity_star = _star_state(gamma)
    similarity = (np.asarray(x) - discontinuity) / time
    density = np.empty_like(similarity)
    velocity = np.empty_like(similarity)
    pressure = np.empty_like(similarity)

    head_speed = -sound_l
    sound_star_l = sound_l * (p_star / p_l) ** ((gamma - 1.0) / (2.0 * gamma))
    tail_speed = velocity_star - sound_star_l
    density_star_l = rho_l * (p_star / p_l) ** (1.0 / gamma)
    shock_speed = sound_r * np.sqrt(
        (gamma + 1.0) / (2.0 * gamma) * p_star / p_r
        + (gamma - 1.0) / (2.0 * gamma)
    )
    density_star_r = rho_r * (
        (p_star / p_r + (gamma - 1.0) / (gamma + 1.0))
        / ((gamma - 1.0) / (gamma + 1.0) * p_star / p_r + 1.0)
    )

    left_data = similarity <= head_speed
    fan = (similarity > head_speed) & (similarity <= tail_speed)
    left_star = (similarity > tail_speed) & (similarity <= velocity_star)
    right_star = (similarity > velocity_star) & (similarity < shock_speed)
    right_data = similarity >= shock_speed

    density[left_data], velocity[left_data], pressure[left_data] = rho_l, 0.0, p_l
    fan_factor = 2.0 / (gamma + 1.0) - (
        (gamma - 1.0) * similarity[fan] / ((gamma + 1.0) * sound_l)
    )
    density[fan] = rho_l * fan_factor ** (2.0 / (gamma - 1.0))
    velocity[fan] = 2.0 * (sound_l + similarity[fan]) / (gamma + 1.0)
    pressure[fan] = p_l * fan_factor ** (2.0 * gamma / (gamma - 1.0))
    density[left_star], velocity[left_star], pressure[left_star] = (
        density_star_l,
        velocity_star,
        p_star,
    )
    density[right_star], velocity[right_star], pressure[right_star] = (
        density_star_r,
        velocity_star,
        p_star,
    )
    density[right_data], velocity[right_data], pressure[right_data] = rho_r, 0.0, p_r
    return np.stack((density, velocity, pressure))
