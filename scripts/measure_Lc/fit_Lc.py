#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import h5py, argparse
import numpy as np
from pylab import figure, close

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
Rperp = fi['Rxx_perp'][...] + fi['Ryy_perp'][...]
Rpara = fi['Rxx_para'][...] + fi['Ryy_para'][...]

Rperp_avr = Rperp.mean(axis=0)
Rperp_std = Rperp.std(axis=0)

Rpara_avr = Rpara.mean(axis=0)
Rpara_std = Rpara.std(axis=0)


fig = figure(1, figsize=(6,4))
ax  = fig.add_subplot(111)

ax.plot(dr, Rpara_avr, '-o', ms=3)
inf, sup = Rpara_avr-Rpara_std, Rpara_avr+Rpara_std
ax.fill_between(dr, inf, sup, facecolor='gray', alpha=0.5)

ax.grid(True)
fig.savefig(pa.figname, dpi=135, bbox_inches='tight')
close(fig)

fi.close()
#EOF
