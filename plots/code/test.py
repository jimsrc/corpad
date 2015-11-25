from pylab import *
from numpy import *

dir_post	= '../../post'
Nm		= 128
perc_slab       = 0.20
sig             = 1e0
Lc_2d           = 1e-2
Lc_slab         = 1e-2
fname_inp	= '%s/diff.pars_Nm%03d_slab%1.2f_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e.dat' % (dir_post, Nm, perc_slab, sig, Lc_2d, Lc_slab)

data		= loadtxt(fname_inp, unpack=True)
EKs		= data[0]
RLs		= data[1]
kperp		= data[2]
kzz		= data[3]
MFPs_perp	= data[4]
MFPs_parall	= data[5]

AUincm  = 1.5e13        		# [cm]
fig	= figure(1, figsize=(4, 7))
ax	= []
ax	+= [fig.add_subplot(311)]
ax	+= [fig.add_subplot(312)]
ax	+= [fig.add_subplot(313)]

adim_RLs	= RLs/(AUincm*Lc_slab)
adim_MFPs_perp	= MFPs_perp/(AUincm*Lc_slab)
adim_MFPs_parall= MFPs_parall/(AUincm*Lc_slab)
ax[0].plot(adim_RLs, adim_MFPs_perp, '-o', c='blue')
ax[1].plot(adim_RLs, adim_MFPs_parall, '-o', c='green')
ax[2].plot(adim_RLs, MFPs_perp/MFPs_parall, '-o', c='black')

TITLE1	= '$L_c^{2D}[AU]: %1.1e$   $L_c^{slab}[AU]: %1.1e$' % (Lc_2d, Lc_slab)
TITLE0	= '$Nm:%d$    $perc\_slab$: $%1.2f$    $(\sigma/Bo)^2$: $%1.1e$' %  (Nm, perc_slab, sig)
TITLE 	= TITLE0 +'\n'+ TITLE1
ax[0].set_title(TITLE)
ax[0].set_ylabel('$\lambda_{\perp} / L_c^{slab}$')
ax[1].set_ylabel('$\lambda_{\parallel} / L_c^{slab}$')
ax[2].set_ylabel('$\lambda_{\perp} / \lambda_{\parallel}$')
ax[2].set_xlabel('$R_L / L_c^{slab}$')
for i in range(3):
	ax[i].set_xscale('log')
	ax[i].set_yscale('log')
	ax[i].grid()

fname_fig	= './test.png'
savefig(fname_fig, format='png', dpi=135, bbox_inches='tight')
close()
