#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import os, argparse
import h5py

# retrieve args
parser = argparse.ArgumentParser(
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-LcS', '--LcSlab', 
type=float, 
default=0.01,  # [AU]
help='integral scale of Slab model'
)
parser.add_argument(
'-xi', '--xi', 
type=float, 
default=1.0,  # [1]
help='ratio Lc2D/LcSlab'
)
parser.add_argument(
'-Bo', '--Bo', 
type=float, 
default=5e-5, # [G] at solar wind, 1AU
help='background mean field Bo (in Gauss)'
)
parser.add_argument(
'-fs', '--fnames', 
type=str, 
nargs='+',
default=None, # [G] at solar wind, 1AU
help='list of .h5 input files'
)
pa = parser.parse_args()

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

for fname_inp in pa.fnames:
    #fname_inp = '../out/h.004/post.h5'
    f = h5py.File(fname_inp, 'r')
    # `psim/Lc_slab` is in units of Larmor radii
    RloLc = 1./f['psim/Lc_slab'].value # [1]
    o = get_phys(RloLc=RloLc, Lc_slab=pa.LcSlab, Bo=pa.Bo)
    print fname_inp, o['Rl']/AUincm, '%e'%o['Ek']



#EOF
