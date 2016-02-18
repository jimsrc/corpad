from pylab import *
from numpy import *

dir_figs	= '../../plots'
dir_post        = '../../../../post'
Nm              = 128 
perc_slab       = 0.20
sig             = 1e0 
Lc_2d           = 1e-2
Lc_slab         = 1e-2
fname_inp       = '%s/diff.pars_Nm%03d_slab.%1.2f_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e.dat' % (dir_post, Nm, perc_slab, sig, Lc_2d, Lc_slab)
fname_theory	= '../../../../../theory/WNLT/composite/mfpParall.vs.Energy/mfp.profiles/xminSlab.1.0e-06_xmaxSlab1.0e+05_xmin2d.1.0e-04_xmax2d.1.0e+03/mfpParallCompos_slab%1.2f_sig%1.2f_Lc2D.%1.1e_LcSlab.%1.1e.dat' % (perc_slab, sig, Lc_2d, Lc_slab)
print " leido: %s" % fname_inp
print " leido: %s" % fname_theory
#------------------------ simulations
data            = loadtxt(fname_inp, unpack=True)
EKs             = data[0]	# [eV]
RLs             = data[1]	# [cm]
kperp           = data[2]	# [cm2/s]
kzz             = data[3]	# [cm2/s]
MFPs_perp       = data[4]	# [cm]
MFPs_parall     = data[5]	# [cm]
#------------------------ theory WNLT
data_theory	= loadtxt(fname_theory, unpack=True)
theory_Rnorm	= data_theory[0] # [1]
theory_mfpParall = data_theory[1] # [AU]
#-------------------------------------

AUincm  = 1.5e13                        # [cm]
fig     = figure(1, figsize=(6, 4))
ax	= fig.add_subplot(111)

adim_RLs        = RLs/(AUincm*Lc_slab)
adim_MFPs_perp  = MFPs_perp/(AUincm*Lc_slab)
adim_MFPs_parall= MFPs_parall/(AUincm*Lc_slab)

ax.plot(adim_RLs, adim_MFPs_parall, '-o', c='black', label='simulations')
ax.plot(theory_Rnorm, theory_mfpParall/Lc_slab, '-', c='blue', label='WNLT approximation')
TITLE1  = '$L_c^{2D}[AU]: %1.1e$   $L_c^{slab}[AU]: %1.1e$' % (Lc_2d, Lc_slab)
TITLE0  = '$Nm:%d$    $\sigma_{slab}^2/\sigma^2$: $%1.2f$    $(\sigma/Bo)^2$: $%1.1e$' %  (Nm, perc_slab, sig)
TITLE   = TITLE0 +'\n'+ TITLE1 + '\n'
ax.set_title(TITLE, fontsize=17.)
ax.set_xlabel('$R_L / L_c^{slab}$', fontsize=19.)
ax.set_ylabel('$\lambda_{\perp} / L_c^{slab}$', fontsize=19.)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(3e-3, 2e0)
ax.legend(loc='upper left')
ax.grid()

fname_fig       = '%s/mfpParall.vs.Rnorm_Lc2d.%1.1e_LcSlab.%1.1e_sig.%1.2f_Nm%03d_slab.%1.2f.png' % (dir_figs, Lc_2d, Lc_slab, sig, Nm, perc_slab)
savefig(fname_fig, format='png', dpi=135, bbox_inches='tight')
print " ----> generamos: %s" % fname_fig
close()
##
