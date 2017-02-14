#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import h5py, argparse
import numpy as np
from pylab import figure, close
import shared.fhist2d as fhist

#--- retrieve args
parser = argparse.ArgumentParser(
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-i', '--input',
type=str,
default='./test.h5',
help='HDF5 input for R-function profiles',
)
parser.add_argument(
'-fig', '--figname',
type=str,
default='./test.png',
help='output figure',
)
pa = parser.parse_args()

fi = h5py.File(pa.input, 'r')
#--- read data
Lc    = fi['psim/Lc_slab'].value  # we assume: slab & 2d with the same Lc
Nrlz  = fi['psim/Nrlz'].value
dr    = fi['dr'][...] # spatial displacements
dth   = fi['dth'][...] # spatial displacements
Rxx   = fi['Rxx'][...]
Ryy   = fi['Ryy'][...]
R     = 0.5*(Rxx + Ryy)

Ravr  = R.mean(axis=0)
Rstd  = R.std(axis=0)

#--- figure
fig = figure(1, figsize=(6,4))
ax  = fig.add_subplot(111)

r_para = np.outer(dr, np.cos(dth)).T
r_perp = np.outer(dr, np.sin(dth)).T

fig, ax = fhist.contour_2d(fig, ax, r_perp, r_para, Ravr.T, hscale='linear', cb_label='Rfunc', vmin=None, vmax=None)

ax.set_xlabel('$r_\perp / L_c$')
ax.set_ylabel('$r_\parallel / L_c$')
fig.savefig(pa.figname, dpi=135, bbox_inches='tight')
close(fig)


#EOF
