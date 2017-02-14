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
Rperp = .5*(fi['Rxx_perp'][...] + fi['Ryy_perp'][...])
Rpara = .5*(fi['Rxx_para'][...] + fi['Ryy_para'][...])

Rperp_avr = Rperp.mean(axis=0)
Rperp_std = Rperp.std(axis=0)

Rpara_avr = Rpara.mean(axis=0)
Rpara_std = Rpara.std(axis=0)


#--- figure
fig = figure(1, figsize=(9,4))

#--- fit parallel
cc = dr<1.3 # select interval for the fit
x = dr[cc]
y = np.log(Rpara_avr[cc])
err = Rpara_std/np.sqrt(Nrlz) # error
m, b = np.polyfit(x, y, deg=1, w=1./err[cc], cov=False)
lc_para = -1./m # fitted correlation length
#--- fig parallel
ax  = fig.add_subplot(1,2,1)
ax.set_title("Nrlz:%d" % fi['psim/Nrlz'].value)
ax.plot(dr, Rpara_avr, '-ok', ms=3)
inf, sup = Rpara_avr-err, Rpara_avr+err
ax.fill_between(dr, inf, sup, facecolor='gray', alpha=0.5)
# plot the fit
ax.plot(x, np.exp(m*x+b), '-r', lw=3, label='$\lambda_c=%1.3g$'%lc_para, alpha=.6)
ax.set_xlabel('$r_\parallel$')
ax.set_ylabel('$R_\perp(r_\parallel$)')
ax.legend(loc='best')
ax.grid(True)
#ax.set_yscale('log')


#--- fit perp
cc = dr<1.3 # select interval for the fit
x = dr[cc]
y = np.log(Rperp_avr[cc])
err = Rperp_std/np.sqrt(Nrlz) # error
m, b = np.polyfit(x, y, deg=1, w=1./err[cc], cov=False)
lc_perp = -1./m # fitted correlation length
#--- fig perp
ax  = fig.add_subplot(1,2,2)
ax.set_title("Nrlz:%d" % fi['psim/Nrlz'].value)
ax.plot(dr, Rperp_avr, '-ok', ms=3)
inf, sup = Rperp_avr-err, Rperp_avr+err
ax.fill_between(dr, inf, sup, facecolor='gray', alpha=0.5)
# plot the fit
ax.plot(x, np.exp(m*x+b), '-r', lw=3, label='$\lambda_c=%1.3g$'%lc_perp, alpha=.6)
ax.set_xlabel('$r_\perp$')
ax.set_ylabel('$R_\perp(r_\perp$)')
ax.legend(loc='best')
ax.grid(True)
#ax.set_yscale('log')

fig.tight_layout(w_pad=1.)
fig.savefig(pa.figname, dpi=135, bbox_inches='tight')
close(fig)

fi.close()
#EOF
