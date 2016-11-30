#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
from pylab import figure, show, close
import h5py, argparse, os

#--- universal constants
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
