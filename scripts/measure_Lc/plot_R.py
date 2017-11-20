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

print(' [+] reading %s\n'%pa.input)
fi = h5py.File(pa.input, 'r')
#--- read data
Lc    = fi['psim/Lc_slab'].value  # we assume: slab & 2d with the same Lc
Nrlz  = fi['psim/Nrlz'].value
xi    = fi['psim/xi'].value
dr    = fi['dr'][...] # spatial displacements
dth   = fi['dth'][...] # spatial displacements
Rxx   = fi['Rxx'][...]
Ryy   = fi['Ryy'][...]
R     = 0.5*(Rxx + Ryy)

print(' [*] averaging over %d realizations...\n' % Nrlz)
Ravr  = R.mean(axis=0)
Rstd  = R.std(axis=0)

#--- figure R vs (r,th)
fig = figure(1, figsize=(11,4))
ax  = fig.add_subplot(1,2,1)

r_para = np.outer(dr, np.cos(dth)).T
r_perp = np.outer(dr, np.sin(dth)).T

print(' [*] plotting R vs (r,th)\n')
fig, ax = fhist.contour_2d(fig, ax, r_perp, r_para, Ravr.T, hscale='linear', cb_label='Rfunc', vmin=None, vmax=None)

ax.set_title("$Nrlz = %d$" % fi['psim/Nrlz'].value)
ax.set_xlabel('$r_\perp / L_c$')
ax.set_ylabel('$r_\parallel / L_c$')
#fig.savefig(pa.figname, dpi=135, bbox_inches='tight')
#close(fig)

#--- figure R vs r
fig = figure(1, figsize=(6,4))
ax  = fig.add_subplot(1,2,2)

print(' [*] averaging over theta...\n')
Ravr_r = (R.mean(axis=2)).mean(axis=0)
Rstd_r = (R.mean(axis=2)).std(axis=0)
err = Rstd_r/np.sqrt(Nrlz) # error

print(' [*] plotting R vs r (with error bands) ...\n')
ax.plot(dr, Ravr_r, '-ok', ms=3, label='avrg over $\\theta$')
inf, sup = Ravr_r-err, Ravr_r+err
ax.fill_between(dr, inf, sup, facecolor='gray', alpha=0.5)

#--- fit 
print(' [*] performing exponential fit...\n')
cc = dr<5. # select interval for the fit
x = dr[cc]
y = np.log(Ravr_r[cc])
err = Rstd_r/np.sqrt(Nrlz) # error
m, b = np.polyfit(x, y, deg=1, w=1./err[cc], cov=False)
lcorr = -1./m # fitted correlation length
ax.plot(x, np.exp(m*x+b), '-r', lw=3, label='$\lambda_c=%1.3g$'%lcorr, alpha=.6)

ax.set_xlabel('$r/L_c$')
ax.set_ylabel('$\langle R \\rangle_\\theta$')
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid(True)
TITLE = '$L_c = %g$;   $\\xi=L_c^{2d}/L_c^{slab} = $%g' % (Lc,xi)
ax.set_title(TITLE)

print(' [*] saving to: %s\n' % pa.figname)
fig.tight_layout(w_pad=1.)
fig.savefig(pa.figname, dpi=135, bbox_inches='tight')
close(fig)
fi.close()

#--- let's append the fit-info
print(' [*] appending fitted parameters to: %s\n' % pa.input)
fi = h5py.File(pa.input, 'r+')
if 'fit' in fi:
    fi.pop('fit') # over-write dataset
fi['fit/lcorr'] = lcorr
fi.close()

#EOF
