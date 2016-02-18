from pylab import *
from numpy import *

dir_figs	= '../../plots'
dir_post        = '../../../../post'
dir_theory	= '../../../../../theory/WNLT/composite'
Nm              = 128 
perc_slab       = 0.20
sig             = 1e0 
Lc_2d           = 3e-3
Lc_slab         = 3e-2
Ek		= 6e5	# [eV]
R		= 5e-3	# [1] radio larmor normalizado con Lc_slab
fname_inp	= '%s/diff.pars_Nm%03d_Ek.%1.1e_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e.dat' % (dir_post, Nm, Ek, sig, Lc_2d, Lc_slab)
fname_theory	= '%s/mfpParall.vs.SlabPerc/mfp.profiles/xminSlab.1.0e-06_xmaxSlab1.0e+05_xmin2d.1.0e-04_xmax2d.1.0e+03/mfp.ParallCompos_R.%1.4f_sig%1.2f_Lc2D.%1.1e_LcSlab.%1.1e.dat' % (dir_theory, R, sig, Lc_2d, Lc_slab)
print " leido: %s" % fname_inp
print " leido: %s" % fname_theory

#------------------------ simulations
data            = loadtxt(fname_inp, unpack=True)
perc_slabs	= data[0]	# [eV]
Kperp           = data[1]	# [cm2/s]
kparall         = data[2]	# [cm2/s]
MFPs_perp       = data[3]	# [cm]
MFPs_parall     = data[4]	# [cm]
#------------------------ theory WNLT
data_theory	= loadtxt(fname_theory, unpack=True)
theory_slabs	= data_theory[0] # [1]
theory_mfpParall= data_theory[1] # [AU]
#-------------------------------------

AUincm  = 1.5e13                        # [cm]
fig     = figure(1, figsize=(6, 4))
ax	= fig.add_subplot(111)

ax.plot(perc_slabs[1:], MFPs_parall[1:]/(AUincm*Lc_slab), '-o', c='black', label='simulations')
ax.plot(theory_slabs, theory_mfpParall/Lc_slab, '-', c='blue', label='WNLT approximation')
TITLE1  = '$L_c^{2D}[AU]: %1.1e$   $L_c^{slab}[AU]: %1.1e$' % (Lc_2d, Lc_slab)
TITLE0  = '$Nm:%d$    $R_L/L_c^{slab}$: $%1.1e$    $(\sigma/Bo)^2$: $%1.1e$' %  (Nm, R, sig)
TITLE   = TITLE0 +'\n'+ TITLE1 + '\n'
ax.set_title(TITLE, fontsize=17.)
ax.set_xlabel('$\sigma_{slab}^2 / \sigma^2$', fontsize=19.)
ax.set_ylabel('$\lambda_{\parallel} / L_c^{slab}$', fontsize=19.)
#ax.set_xscale('log'); 
#ax.set_yscale('log')
ax.set_ylim(0., 5)
ax.legend(loc='upper left')
ax.grid()

fname_fig       = '%s/mfpParall.vs.Slabs_Lc2d.%1.1e_LcSlab.%1.1e_sig.%1.2f_Nm%03d_sig.%1.1e_R.%1.1e.png' % (dir_figs, Lc_2d, Lc_slab, sig, Nm, sig, R)
savefig(fname_fig, format='png', dpi=135, bbox_inches='tight')
print " ----> generamos: %s" % fname_fig
close()
##"""
