"""Checks of the model.  Run:  python -m pytest -q"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from model import atom
from model.ring import (HBAR, Params, coherence, current, evolve, evolve_expm, ground_fraction,
                        initial_state, operators, physicality, record_statistics, site_probabilities)


def run(params, kind="circulating", t=12.0):
    ops = operators(params)
    rhos = evolve(ops, initial_state(kind), np.array([0.0, t]))
    return ops, rhos


def test_hamiltonian_is_hermitian_and_has_81_states():
    H = operators()["H"]
    assert H.shape == (81, 81)
    assert np.allclose(H, H.conj().T)


def test_state_stays_physical():
    ops, rhos = run(Params())
    trace_err, herm_err, min_eig = physicality(rhos[-1])
    assert trace_err < 1e-10 and herm_err < 1e-8 and min_eig > -1e-8
    assert np.allclose(site_probabilities(rhos[-1]), 1 / 3, atol=1e-10)


def test_without_mixing_nothing_changes():
    ops, rhos = run(Params(g=0.0))
    assert coherence(rhos[-1]) == pytest.approx(1.0, abs=1e-8)
    assert current(ops, rhos[-1]) / current(ops, rhos[0]) == pytest.approx(1.0, abs=1e-8)
    assert ground_fraction(ops, rhos[-1]) == pytest.approx(0.0, abs=1e-12)


def test_free_hopping_matches_exact_formula():
    ops = operators(Params(g=0.0))
    times = np.linspace(0, 12, 25)
    rhos = evolve(ops, initial_state("localized"), times)
    p0 = np.array([np.real(np.trace(ops["n"][0] @ r)) for r in rhos])
    assert np.allclose(p0, (5 + 4 * np.cos(3 * times / HBAR)) / 9, atol=1e-8)


def test_two_independent_solvers_agree():
    ops, rhos = run(Params())
    assert np.abs(evolve_expm(ops, rhos[0], 12.0) - rhos[-1]).max() < 1e-9


def test_1s_energy_cancels():
    ops, rhos = run(Params())
    ops2, rhos2 = run(Params(e1s=-37.0))
    assert coherence(rhos[-1]) == pytest.approx(coherence(rhos2[-1]), abs=1e-9)


def test_reference_numbers_at_12_ns():
    ops, rhos = run(Params())
    assert current(ops, rhos[-1]) / current(ops, rhos[0]) == pytest.approx(0.2171, abs=5e-4)
    assert coherence(rhos[-1]) == pytest.approx(0.4110, abs=5e-4)
    assert ground_fraction(ops, rhos[-1]) == pytest.approx(0.3859, abs=5e-4)
    info, joint = record_statistics(rhos[-1])
    assert info == pytest.approx(0.00823, abs=5e-5)
    assert 1 - joint[:, 0].sum() == pytest.approx(0.7070, abs=5e-4)


def test_emission_off_keeps_all_atoms_excited():
    ops, rhos = run(Params(gamma=0.0))
    assert ground_fraction(ops, rhos[-1]) == pytest.approx(0.0, abs=1e-12)
    assert current(ops, rhos[-1]) / current(ops, rhos[0]) == pytest.approx(0.4363, abs=5e-4)


def test_atoms_run_out():
    ops, rhos = run(Params(), t=300.0)
    assert ground_fraction(ops, rhos[-1]) > 0.999


def test_dipole_matrix_element():
    assert atom.dipole_length_orbital() == pytest.approx(-3.0, abs=1e-9)
    assert atom.dipole_length_fine_structure() == pytest.approx(np.sqrt(3), abs=1e-9)


def test_quench_rate_weak_field_and_ceiling():
    V = atom.mixing_energy(1e3)
    assert V == pytest.approx(0.0917, abs=1e-4)
    assert atom.quench_rate_weak_field(V) * 1e9 == pytest.approx(2.744e5, rel=1e-3)
    assert atom.quench_rate_exact(V) == pytest.approx(atom.quench_rate_weak_field(V), rel=5e-3)
    assert atom.quench_rate_exact(1e4) == pytest.approx(atom.GAMMA_2P / 2, rel=1e-3)
    for field in np.logspace(1, 9, 30):
        assert atom.quench_rate_exact(atom.mixing_energy(field)) <= atom.GAMMA_2P / 2 * (1 + 1e-9)
