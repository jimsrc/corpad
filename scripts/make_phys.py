#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import os, argparse
import h5py
import funcs as ff
#from src_Bmodel import Bmodel # turb model
from shared.funcs import calc_Rlarmor, Bo_parker, Lc_memilia
from pylab import figure, show, close

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

ro      = 1.0  # [AU]
Bo      = Bo_parker(r=ro)  # [G]
Lc_slab = Lc_memilia(r=ro) # [AU]
#Rl = calc_Rlarmor(
#    rigidity=1e9,       # [V] rigidity
#    Bo=Bo,              # [G] Bo-field
#    )/AUincm            # [AU]
print " ro: %e" % ro 
print " Bo: %e" % Bo
print " Lc: %e" % Lc_slab 
print " -------------------- "

Ek, lparall, lperp = [], [], []
RloLc = []
for fname_inp in pa.fnames:
    #fname_inp = '../out/h.004/post.h5'
    f = h5py.File(fname_inp, 'r')
    # `psim/Lc_slab` is in units of Larmor radii
    _RloLc = 1./f['psim/Lc_slab'].value # [1]
    o = ff.get_phys(RloLc=_RloLc, Lc_slab=Lc_slab, Bo=Bo)
    print fname_inp, o['Rl']/AUincm, '%e'%o['Ek']
    #--- set B-turbulence model
    #Rl = o['Rl']/AUincm  # [AU]
    RloLc   += [ _RloLc ]
    Ek      += [ o['Ek'] ]
    lparall += [ f['pfit/lparall'].value ]
    lperp   += [ f['pfit/lperp'].value ]

#import pdb; pdb.set_trace()
fig = figure(1, figsize=(6,4))
ax  = fig.add_subplot(111)
ax2 = ax.twinx()
ax.plot(Ek, lparall, '-ob', label='parall')
ax2.plot(Ek, lperp, '-or', label='perp')
ax.set_yscale('log')
ax2.set_yscale('log')
ax.set_xscale('log')
ax.grid()
ax.legend(handles=ax.get_legend_handles_labels()[0], loc='best')
fig.savefig('./test.png', dpi=138, bbox_inches='tight')
close(fig)

#EOF
