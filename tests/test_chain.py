"""Checks of the long-chain code.  Run:  python -m pytest -q"""
import os
import sys

import numpy as np
import pytest
from scipy.special import jv

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from model.chain import (HBAR, ChainParams, SectorCache, exact_density_matrix, free_hopping,
                         no_emission, trajectory)


def average(p, start, times, n, seed0=0):
    cache = SectorCache(p)
    runs = np.array([trajectory(p, start, times, np.random.default_rng(seed0 + k), cache)[0] for k in range(n)])
    return runs.mean(axis=0), runs.std(axis=0) / np.sqrt(n)


def test_free_hopping_is_a_bessel_function():
    N, A, times = 41, 0.1, np.linspace(0, 40, 9)
    p = free_hopping(N, A, N // 2, times)
    x = np.arange(N) - N // 2
    assert np.allclose(p, jv(x[None, :], 2 * A * times[:, None] / HBAR) ** 2, atol=1e-6)


def test_emission_off_matches_exact_solution():
    p = ChainParams(N=3, A=0.3, gamma=0.0, max_2p=3)
    times = np.linspace(0, 20, 6)
    assert np.allclose(no_emission(p, 1, times), exact_density_matrix(p, 1, times), atol=1e-6)


def test_trajectories_match_exact_lindblad_solution():
    p = ChainParams(N=3, A=0.3, max_2p=3)
    times = np.linspace(0, 30, 7)
    mean, err = average(p, 1, times, 400)
    exact = exact_density_matrix(p, 1, times)
    assert np.all(np.abs(mean - exact) < 4 * err + 0.01)


def test_two_simultaneous_2p_atoms_are_enough():
    """Evolution without a photon: at most two 2p atoms gives the same state as at most four."""
    out = []
    for m in (2, 4):
        p = ChainParams(N=7, max_2p=m)
        sector = SectorCache(p)[tuple(range(7))]
        psi = np.zeros(sector.dim, dtype=complex)
        psi[3 * sector.nsub] = 1.0
        out.append(sector.site_probabilities(sector.propagate(psi, 20.0)))
    assert np.abs(out[0] - out[1]).max() < 1e-3


def test_level_repulsion_shift():
    p = ChainParams()
    shift = (-p.delta + np.sqrt(p.delta ** 2 + 4 * p.g ** 2)) / 2
    assert shift == pytest.approx(0.523, abs=1e-3)
    assert np.linalg.eigvalsh(np.array([[0.0, p.g], [p.g, -p.delta]])).max() == pytest.approx(shift)


def test_electron_stays_near_its_start_with_slow_hopping():
    p = ChainParams(N=11, A=0.1, max_2p=2)
    times = np.array([0.0, 60.0])
    mean, _ = average(p, 5, times, 40)
    kept = mean[-1, 3:8].sum()
    free = sum(jv(n, 2 * p.A * 60 / HBAR) ** 2 for n in range(-2, 3))
    assert kept > 0.6 and free < 0.2
