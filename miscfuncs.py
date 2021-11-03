#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import scipy as sp
import qutip as qt
import numba as nb
from math import erf, lgamma
import ctypes

from numba.extending import get_cython_function_address
addr = get_cython_function_address("scipy.special.cython_special",
                                   "eval_hermitenorm")
functype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_int, ctypes.c_double)
eval_hermitenorm = functype(addr)

def displacement_fock(alpha, dim):
    def dmn(m, n, alpha):
        logfac = (0.5 * (lgamma(n + 1) - lgamma(m + 1))
                  - 0.5 * np.abs(alpha) ** 2 + (m - n) * np.log(np.abs(alpha)))
        if logfac < -700:
            return 0.0
        return (np.exp(logfac) * np.exp(1j*np.angle(alpha) * (m - n))
                * sp.special.eval_genlaguerre(n, m - n, np.abs(alpha) ** 2))
    
    if type(dim) is tuple:
        (m, n) = dim
        if m > n:
            return dmn(m, n, alpha)
        else:
            return dmn(n, m, np.conj(-alpha))
    
    elif type(dim) is int:
        if alpha == 0:
            return qt.qeye(dim)
        else:   
            dmat = np.zeros((dim, dim), dtype=np.complex_) 
            for m in range(0, dim):
                for n in range(0, dim):
                    if m > n:
                        d = dmn(m, n, alpha)
                    else:
                        d = dmn(n, m, np.conj(-alpha))
                    dmat[m, n] = d
            return qt.Qobj(dmat)
        
@nb.njit(fastmath=True)
def hermint_4_numba(d, a, phase):
    i1 = np.ceil(np.sqrt(2 * d + 1) / a) + 1
    imax = int(i1 + np.mod(i1, 2))
    M = np.zeros((d, d), dtype=np.complex_)

    for s in range(imax):
        M[0, 0] += ((-1) ** s) * erf(a * (s + 0.5))
    M[0, 0] += erf(imax * a) / 2
    M[0, 0] *= 2

    for m in range(2, d, 2):
        sum = 0.0
        logfacm = np.log(2) + m / 2 * np.log(2) - lgamma(m + 1) / 2
        for s in range(imax + 1):
            sa = (s + 1 / 2) * a
            x = np.exp(- sa ** 2) * 2.0 ** (-(m - 1) / 2)                 * eval_hermitenorm(m - 1, np.sqrt(2) * sa)
            sum += ((-1) ** s) * np.exp(logfacm) * x
        sum *= (-np.exp(1j * m * phase ) / np.sqrt(np.pi))
        M[m, 0] = sum
        M[0, m] = np.conj(M[m, 0])

    for n in range(1, d):
        for m in range(np.mod(n + 1, 2) + 1, n + 1, 2):
            sum = 0.0
            logfac = np.log(2) + m / 2 * np.log(2) - lgamma(m + 1) / 2                      + n / 2 * np.log(2) - lgamma(n + 1) / 2
            for s in range(imax + 1):
                sa = (s + 1 / 2) * a
                x = np.exp(- sa ** 2 / 2) * 2.0 ** (-m / 2)                     * eval_hermitenorm(m, np.sqrt(2) * sa)
                y = np.exp(- sa ** 2 / 2) * 2.0 ** (-(n - 1) / 2)                     * eval_hermitenorm(n - 1, np.sqrt(2) * sa)
                sum += ((-1) ** s) * np.exp(logfac) * x * y / np.sqrt(np.pi)
            sum *= -np.exp(1j * (m - n) * phase)
            add = np.sqrt(m / n) * M[m - 1, n - 1]
            M[m, n] = sum + add
            z = np.abs(M[m, n])
            M[n, m] = np.conj(M[m, n])
    return M

def gkp_measure_ops(dim, alpha, beta):
    gamma = alpha + beta
    x_measure = hermint_4_numba(dim, np.pi/(np.abs(alpha)*np.sqrt(2)),
                            np.angle(alpha) + np.pi/2)
    y_measure = hermint_4_numba(dim, np.pi/(np.abs(gamma)*np.sqrt(2)),
                            np.angle(gamma) - np.pi/2)
    z_measure = hermint_4_numba(dim, np.pi/(np.abs(beta)*np.sqrt(2)),
                            np.angle(beta) - np.pi/2)
    return qt.Qobj(x_measure), qt.Qobj(y_measure), qt.Qobj(z_measure)

def fock_to_pos(psi, q, phi=0.0):
    if isinstance(psi, qt.Qobj):
        psin = psi.data.toarray().flatten()
    else:
        psin = psi
    psin = np.asarray(psin)
    if np.isscalar(q):
        qv = np.array([q])
    else:
        qv = np.asarray(q)
    return _fock_to_posv(psin, qv, np.array([phi]))


@nb.guvectorize([(nb.complex128[:], nb.float64[:], nb.float64[:],
                  nb.complex128[:])],
                '(n),(),()->()', nopython=True)
def _fock_to_posv(psin, qv, phiv, psiqv):
    def c(m, q):
        logcm = - 0.5 * lgamma(m + 1)
        return np.exp(logcm) * eval_hermitenorm(m, np.sqrt(2) * q)

    q = qv[0]
    phi = phiv[0]
    cm2 = c(0, q) # c0
    cm1 = c(1, q) # c1
    psiq = cm2 * psin[0] + cm1 * np.exp(-1j*phi) * psin[1]
    for m in range(2, len(psin)):
        cm = np.sqrt(2) * q/np.sqrt(m) * cm1 - np.sqrt((m-1)/m) * cm2
        cm2 = cm1
        cm1 = cm
        psiq += cm * np.exp(-1j*phi*m) * psin[m]
    psiq = psiq * np.exp(-q ** 2 / 2) / (np.pi) ** (1 / 4)
    psiqv[0] = psiq

