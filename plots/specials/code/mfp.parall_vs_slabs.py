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

fig	= figure(1, figsize=(6, 4))
ax	= fig.add_subplot(111)

ax.plot(slab_percs, mfp_parall/AUincm, '-o', c='red', label='our simulations')

ax.grid()
ax.set_xlabel('$\sigma_{slab}^2 / \sigma^2$', fontsize=17)
ax.set_ylabel('$\lambda_\parallel$ [AU]', fontsize=17)
ax.set_ylim(0., 0.25)
ax.legend()

fname_fig = '../mfp.parall_vs_slabs.png'
savefig(fname_fig, dpi=100, format='png', bbox_inches='tight')
close()
