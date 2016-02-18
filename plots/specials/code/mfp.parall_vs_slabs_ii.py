#!/usr/bin/env ipython
from pylab import *
from numpy import *

AUincm	= 1.5e13		# [cm]
fname	= '../../../post/diff.pars_Nm128_Ek.6.0e+05_sig.1.0e+00_Lc2d.3.0e-03_LcSlab.3.0e-02.dat'
data 	= loadtxt(fname, unpack=True)

slab_percs	= data[0][1:]
kperp		= data[1][1:]
kparall		= data[2][1:]
mfp_perp	= data[3][1:]
mfp_parall	= data[4][1:]
RL          = 1.494646E-04*AUincm # [cm] ;Energ K [eV]= 600000; Bo=5nT; OMEGA [s^-1]= 4.783103E-01  

fig	= figure(1, figsize=(6, 4))
ax	= fig.add_subplot(111)

#ax.plot(slab_percs, mfp_parall/AUincm, '-o', c='red', label='our simulations')
ax.plot(mfp_parall/RL, mfp_perp/mfp_parall , '-o', c='red', label='our simulations')

ax.grid()
ax.set_xlabel('$\lambda_\parallel / R_L$', fontsize=17)
ax.set_ylabel('$\lambda_\perp / \lambda_\parallel$', fontsize=17)
#ax.set_ylim(0., 0.25)
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()

#fname_fig = '../mfp.parall_vs_slabs.png'
fname_fig = './test.png'
savefig(fname_fig, dpi=100, format='png', bbox_inches='tight')
close()
