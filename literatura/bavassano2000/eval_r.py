#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plots sigmaBoB :D
"""
from pylab import *
import numpy as np
from scipy.interpolate import (
    splrep,  # spline interpolation
    splev    # spline evaluator
)
from gg import r as r_eb
from gg import eb 
from Bparker import Bparker as BP

fnm_lam = './fig01_bavassano_etal00_lambda.dat'
fnm_rad = './fig01_bavassano_etal00_radius.dat'

tl, l = np.loadtxt(fnm_lam).T
tr, r = np.loadtxt(fnm_rad).T


time = np.linspace(95.2, 96.6, 500)
#--- interp lambda
tck = splrep(tl, l, s=0)
lam = splev(time, tck, der=0) # der: order of derivative to compute (less than "k")
#--- interp radius
tck = splrep(tr, r, s=0)
rad = splev(time, tck, der=0)
#--- interp lambda=lambda(r)
tck  = splrep(rad, lam, s=0)
llam = splev(r_eb, tck, der=0)
th   = 90.0 - llam # [deg]

th[:2] = np.nan
p = np.zeros(3, dtype=np.float32)
B, Bm = [], []
for p[0], p[1] in zip(r_eb, th):
    br, bt, bp = BP.return_B(p)
    bmod = np.sqrt(br*br+bt*bt+bp*bp)
    B += [[ br, bt, bp, bmod ]]

B  = np.array(B)
#rho = *(1.0/r_eb*r_eb)
#Va = B[:,-1]/rho          # Parker-Alfven veloc.
Va = (40.0)*(B[:,-1]/(5e-5))*r_eb   # [km/s] Parker-Alfven veloc.
sigmaBoB = eb/(Va*Va)

if __name__=='__main__':
    #--- figs
    fig = figure(1, figsize=(6,4))
    ax  = fig.add_subplot(111)
    #ax2 = ax.twinx()

    #--- eb
    ax.plot(r_eb, sigmaBoB, '-*', label='$e_b/V_A^2$')
    #ax.set_xscale('log')
    ax.set_xlim(0.3, 5.)
    ax.grid()

    ax.legend(loc='lower left')
    savefig('./eb_o_Bparker.png', dpi=300, bbox_inches='tight')
    close()
#EOF
