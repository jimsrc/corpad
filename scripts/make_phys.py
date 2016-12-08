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
'-xi', '--xi', 
type=float, 
default=1.0,  # [1]
help='ratio Lc2D/LcSlab'
)
parser.add_argument(
'-fs', '--fnames', 
type=str, 
nargs='+',
default=None, # [G] at solar wind, 1AU
help='list of .h5 input files'
)
parser.add_argument(
'-ro', '--ro', 
type=float, 
nargs='+',
default=None, # [G] at solar wind, 1AU
help="""
list of heliodistances r (in AU), corresponding \
with the input filenames.
"""
)
parser.add_argument(
'-fig', '--figname', 
type=str, 
default='./test.png',  # [1]
help='output figure filename',
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

lparall = {'adim': [], 'phys': []}
lperp   = {'adim': [], 'phys': []}
RloLc = []
r = []
for fname_inp, ro in zip(pa.fnames, pa.ro):
    Bo      = Bo_parker(r=ro)  # [G]
    Lc_slab = Lc_memilia(r=ro) # [AU]
    #Rl = calc_Rlarmor(
    #    rigidity=1e9,       # [V] rigidity
    #    Bo=Bo,              # [G] Bo-field
    #    )/AUincm            # [AU]

    f = h5py.File(fname_inp, 'r')
    # `psim/Lc_slab` is in units of Larmor radii
    _RloLc = 1./f['psim/Lc_slab'].value # [1]
    o = ff.get_phys(RloLc=_RloLc, Lc_slab=Lc_slab, Bo=Bo)
    print fname_inp, o['Rl']/AUincm, '%e'%o['Ek']
    #--- set B-turbulence model
    #Rl = o['Rl']/AUincm  # [AU]
    RloLc   += [ _RloLc ]
    lparall['adim'] += [ f['pfit/lparall'].value ]
    lparall['phys'] += [ f['pfit/lparall'].value*Lc_slab ]
    lperp['adim']   += [ f['pfit/lperp'].value ]
    lperp['phys']   += [ f['pfit/lperp'].value*Lc_slab ]
    r       += [ ro ]

# convert to np.array
for nm in lparall.keys():
    lparall[nm] = np.array(lparall[nm])
for nm in lperp.keys():
    lperp[nm] = np.array(lperp[nm])

#import pdb; pdb.set_trace()
fig = figure(1, figsize=(6,4))
ax  = fig.add_subplot(111)
ax2 = ax.twinx()

ax.plot(r, lparall['phys'], '-ob', label='parall')
ax2.plot(r, lperp['phys'], '-or', label='perp')

#ax.plot(r, RloLc, '-ob', label='parall')

#ax.plot(r, lperp['adim']/lparall['adim'], '-ob', label='parall')

#ax.set_yscale('log')
#ax2.set_yscale('log')
#ax.set_xscale('log')
ax.grid()
ax.legend(
    handles= ax.get_legend_handles_labels()[0] +\
            ax2.get_legend_handles_labels()[0],
    loc='best'
)
fig.savefig(pa.figname, dpi=138, bbox_inches='tight')
print " --> we generated: "+pa.figname
close(fig)

#--- generate ASCII
# list of Lc's
lc = [Lc_memilia(r=ro) for ro in r] # [AU]
rl = [calc_Rlarmor(rigidity=1.69604E+09, Bo=Bo_parker(ro))/AUincm for ro in r] # [AU]
data_o = np.array([r, lparall['phys'], lperp['phys'], lc, rl]).T
fname_out = pa.figname.replace('png','txt')
np.savetxt(fname_out, data_o)

#EOF
