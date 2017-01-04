#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
from pylab import figure, show, close
import h5py, argparse, os
from numpy import power as pow

#--- universal constants
M_PI    = np.pi
M_E     = np.e
c       = 3.0*1e10       # [cm/s] speed of light
q       = 4.8032*1e-10   # [statC] proton charge
nT_in_G = 1.0*1e-5       # [1G=1e5nT]
Eo      = 938272013.0    # [eV] rest mass of proton
AUincm  = 1.5e13         # [cm] 1AU
mo  = (1.6726*1e-24)    # [gr] proton mass
#-----------------------

def get_phys(RloLc, Lc_slab, Bo):
    #---- calculo de la energia
    Rl  = AUincm*Lc_slab*RloLc # [cm]
    """
     note that:
     p*c = Rl*q*Bo
     but,
     p   = gamma*v*mo
     p*c = gamma*beta * mo*c^2
     then,
     gamma*beta = (p*c)/(mo*c^2)
    """
    gb = (Rl*q*Bo)/(mo*c**2)  # [1] gamma*beta
    pc = gb*Eo   # [eV] rest mass of proton
    E  = np.sqrt(Eo**2 + pc**2)  # [eV]
    gamma = E/Eo       # [1]
    Ek = Eo*(gamma-1.)  # [eV]

    out = {
    'Rl'    : Rl, 
    'Ek'    : Ek,
    'gamma' : gamma,
    }
    return out

def spectra2D_GiacaloneAndJokipii(lmin=5e-3, lmax=1e2, Nm=128):
    """
    `k, kmin, kmax`: are assumed in units of Lc (correlation scale)
    """
    kmin    = 2.*M_PI/lmax
    kmax    = 2.*M_PI/lmin
    spec    = np.zeros(Nm, dtype=np.float64)
    k,dk,dV = np.zeros((3,Nm), dtype=np.float64)    # Fourier modes
    g2D     = 8./3.

    # from "defs_turb.cc":
    for i in range(Nm):
        k[i]     = kmin * pow(kmax/kmin, 1.*i/(Nm-1.))
        dk[i]    = k[i] *(pow(kmax/kmin, 1./(Nm-1)) - 1.)
        dV[i]    = 2.*M_PI*k[i]*dk[i]            # volume element in Fourier space
        spec[i]  = 1.0 / (1. + pow(k[i], g2D))
    
    return k, dV, spec


def spectra2D_Shalchi(lmin=5e-3, lmax=1e2, Nm=128):
    """
    `k, kmin, kmax`: are assumed in units of Lc (correlation scale)
    """
    kmin    = 2.*M_PI/lmax
    kmax    = 2.*M_PI/lmin
    spec    = np.zeros(Nm, dtype=np.float64)    # spectra with Area=1.0
    k,dk,dV = np.zeros((3,Nm), dtype=np.float64)    # Fourier modes
    nu      = 5./6.

    # from "defs_turb.cc":
    for i in range(Nm):
        k[i]     = kmin * pow(kmax/kmin, 1.*i/(Nm-1.))
        dk[i]    = k[i] *(pow(kmax/kmin, 1./(Nm-1)) - 1.)
        dV[i]    = 2.*M_PI*k[i]*dk[i]            # volume element in Fourier space
        spec[i]  = 1.0 / (k[i]*pow(1. + k[i]**2, nu))

    return k, dV, spec

#EOF
