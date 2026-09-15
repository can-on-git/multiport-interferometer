
# Portions of this file are derived from Strawberry Fields
# Copyright 2019 Xanadu Quantum Technologies Inc.
#
# Modifications Copyright 2026 Carl Abi Nakad
#
# Modified from:
# https://github.com/XanaduAI/strawberryfields
#
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import numpy as np
from collections import defaultdict

def _rectangular_compact_init(
    U, rtol=1e-12, atol=1e-12
):  # pylint: disable=too-many-statements, too-many-branches
    r"""Rectangular decomposition of a unitary with sMZIs and phase-shifters, as given in FIG. 3 and "The Clements Scheme" section of (arXiv:2104.0756).

    Args:
        U (array): unitary matrix

    Returns:
        dict: A dictionary containing the following items:

        * ``m``: the length of the matrix
        * ``phi_ins``: parameter of the phase-shifter at the beginning of the mode
        * ``sigmas``: parameter of the sMZI :math:`\frac{(\theta_1+\theta_2)}{2}`, where :math:`\theta_{1,2}` are the values of the two internal phase-shifts of sMZI
        * ``deltas``: parameter of the sMZI :math:`\frac{(\theta_1-\theta_2)}{2}`, where :math:`\theta_{1,2}` are the values of the two internal phase-shifts of sMZI
        * ``zetas``: parameter of the phase-shifter at the middle of the mode
        * ``phi_outs``: parameter of the phase-shifter at the end of the mode

    """

    V = U.conj()
    m = U.shape[0]

    phases = {}
    phases["m"] = m
    phases["phi_ins"] = {}  # mode : phi
    phases["deltas"] = {}  # (mode, layer) : delta
    phases["sigmas"] = {}  # (mode, layer) : sigma
    phases["zetas"] = {}  # mode : zeta
    phases["phi_outs"] = {}  # mode : phi

    for j in range(m - 1):
        if j % 2 == 0:
            x = m - 1
            y = j
            phi_j = np.angle(V[x, y + 1]) - np.angle(V[x, y])  # reversed order from paper
            V = V @ P(j, phi_j, m)
            phases["phi_ins"][j] = phi_j
            for k in range(j + 1):
                if V[x, y] == 0:
                    delta = 0.5 * np.pi
                else:
                    delta = np.arctan2(-abs(V[x, y + 1]), abs(V[x, y]))
                n = j - k
                V_temp = V @ M(n, 0, delta, m)
                sigma = np.angle(V_temp[x - 1, y - 1]) - np.angle(V_temp[x - 1, y])
                V = V @ M(n, sigma, delta, m)
                phases["deltas"][n, k] = delta
                phases["sigmas"][n, k] = sigma
                x -= 1
                y -= 1
        else:
            x = m - j - 1
            y = 0
            phi_j = np.angle(V[x - 1, y]) - np.angle(V[x, y])
            V = P(x, phi_j, m) @ V
            phases["phi_outs"][x] = phi_j
            for k in range(j + 1):
                if V[x, y] == 0.0:
                    delta = 0.5 * np.pi
                else:
                    delta = np.arctan2(abs(V[x - 1, y]), abs(V[x, y]))
                V_temp = M(x - 1, 0, delta, m) @ V
                n = m + k - j - 2
                if j != k:
                    sigma = np.angle(V_temp[x + 1, y + 1]) - np.angle(V_temp[x, y + 1])
                else:
                    sigma = 0
                phases["deltas"][n, m - k - 1] = delta
                phases["sigmas"][n, m - k - 1] = sigma
                V = M(n, sigma, delta, m) @ V
                x += 1
                y += 1

    # these next two lines are just to remove a global phase
    zeta = -np.angle(V[0, 0])
    V = V @ P(0, zeta, m)
    phases["zetas"][0] = zeta

    for j in range(1, m):
        zeta = np.angle(V[0, 0]) - np.angle(V[j, j])
        V = V @ P(j, zeta, m)
        phases["zetas"][j] = zeta

    assert np.allclose(V, np.eye(m), rtol=rtol, atol=atol), "decomposition failed"

    return phases

def _absorb_zeta(phases):
    r"""Adjust rectangular decomposition to relocate residual phase-shifters of interferometer to edge-shifters, as given in FIG. 4 and "Relocating residual phase-shifts" section of (arXiv:2104.0756).

    Args:
        phases (dict): output of _rectangular_compact_init

    Returns:
        dict: A dictionary containing the following items:

        * ``m``: the length of the matrix
        * ``phi_ins``: parameter of the phase-shifter at the beginning of the mode
        * ``sigmas``: parameter of the sMZI :math:`\frac{(\theta_1+\theta_2)}{2}`, where :math:`\theta_{1,2}` are the values of the two internal phase-shifts of sMZI
        * ``deltas``: parameter of the sMZI :math:`\frac{(\theta_1-\theta_2)}{2}`, where :math:`\theta_{1,2}` are the values of the two internal phase-shifts of sMZI
        * ``phi_edges``: parameters of the edge phase shifters
        * ``phi_outs``: parameter of the phase-shifter at the end of the mode

    """
    m = phases["m"]
    new_phases = phases.copy()
    del new_phases["zetas"]
    new_phases["phi_edges"] = defaultdict(float)  # (mode, layer) : phi

    if m % 2 == 0:
        new_phases["phi_outs"][0] = phases["zetas"][0]
        for j in range(1, m):
            zeta = phases["zetas"][j]
            layer = m - j
            for mode in range(j, m - 1, 2):
                new_phases["sigmas"][mode, layer] += zeta
            for mode in range(j + 1, m - 1, 2):
                new_phases["sigmas"][mode, layer - 1] -= zeta
            if layer % 2 == 1:
                new_phases["phi_edges"][m - 1, layer] += zeta
            else:
                new_phases["phi_edges"][m - 1, layer - 1] -= zeta
    else:
        for j in range(m):
            zeta = phases["zetas"][j]
            layer = m - j - 1
            for mode in range(j, m - 1, 2):
                new_phases["sigmas"][mode, layer] += zeta
            for mode in range(j + 1, m - 1, 2):
                new_phases["sigmas"][mode, layer - 1] -= zeta
            if layer % 2 == 0:
                new_phases["phi_edges"][m - 1, layer] += zeta
            else:
                new_phases["phi_edges"][m - 1, layer - 1] -= zeta
    return new_phases


def decompose_bell(U, rtol=1e-12, atol=1e-12):
    r"""Rectangular decomposition of a unitary with sMZIs and phase-shifters, as given in FIG. 3+4 and "The Clements Scheme" section of (arXiv:2104.0756).

    Args:
        U (array): unitary matrix

    Returns:
        dict: A dictionary containing the following items:

        * ``m``: the length of the matrix
        * ``phi_ins``: parameter of the phase-shifter at the beginning of the mode
        * ``sigmas``: parameter of the sMZI :math:`\frac{(\theta_1+\theta_2)}{2}`, where :math:`\theta_{1,2}` are the values of the two internal phase-shifts of sMZI
        * ``deltas``: parameter of the sMZI :math:`\frac{(\theta_1-\theta_2)}{2}`, where :math:`\theta_{1,2}` are the values of the two internal phase-shifts of sMZI
        * ``phi_edges``: parameters of the edge phase shifters
        * ``phi_outs``: parameter of the phase-shifter at the end of the mode
    """

    if not U.shape[0] == U.shape[1]:
        raise ValueError("Matrix is not square")

    if not np.allclose(U @ U.conj().T, np.eye(U.shape[0]), rtol=rtol, atol=atol):
        raise ValueError("The input matrix is not unitary")

    phases_temp = _rectangular_compact_init(U, rtol=rtol, atol=atol)
    return _absorb_zeta(phases_temp)

def M(n, sigma, delta, m):
    r"""The symmetric Mach Zehnder interferometer matrix. (Eq 1 of the paper (arXiv:2104.0756).)

    Args:
        n (int): the starting mode of sMZI
        sigma (complex): parameter of the sMZI :math:`\frac{(\theta_1+\theta_2)}{2}`, where :math:`\theta_{1,2}` are the values of the two internal phase-shifts of sMZI
        delta (complex): parameter of the sMZI :math:`\frac{(\theta_1-\theta_2)}{2}`, where :math:`\theta_{1,2}` are the values of the two internal phase-shifts of sMZI
        m (int): the length of the unitary matrix to be decomposed

    Returns:
        array[complex,complex]: the sMZI matrix between n-th and (n+1)-th mode
    """
    mat = np.identity(m, dtype=np.complex128)
    mat[n, n] = np.exp(1j * sigma) * np.sin(delta)
    mat[n, n + 1] = np.exp(1j * sigma) * np.cos(delta)
    mat[n + 1, n] = np.exp(1j * sigma) * np.cos(delta)
    mat[n + 1, n + 1] = -np.exp(1j * sigma) * np.sin(delta)
    return mat


def P(j, phi, m):
    r"""The phase shifter matrix. (Eq 2 of the paper (arXiv:2104.0756).)

    Args:
        j (int): the starting mode of phase-shifter
        phi (complex): parameter of the phase-shifter
        m (int): the length of the unitary matrix to be decomposed

    Returns:
        array[complex,complex]: the phase-shifter matrix on the j-th mode
    """
    mat = np.identity(m, dtype=np.complex128)
    mat[j, j] = np.exp(1j * phi)
    return mat
