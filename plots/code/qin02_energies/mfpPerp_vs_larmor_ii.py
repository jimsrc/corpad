#!/usr/bin/env ipython
from pylab import *
from numpy import *

dir_figs	= '../../mfp.perp'
dir_post        = '../../../post'
dir_theory	= '../../../../theory/wnlt/composite/mfpPerp.vs.Energy/mfp.profiles'
Nm              = 128 
perc_slab       = 0.20
sig             = 1e0 
Lc_2d           = 3e-3
Lc_slab         = 3e-2
fname_inp       = '%s/diff.pars_Nm%03d_slab.%1.2f_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e.dat' % (dir_post, Nm, perc_slab, sig, Lc_2d, Lc_slab)
fname_theory	= '%s/xminSlab.1.0e-06_xmaxSlab1.0e+05_xmin2d.1.0e-04_xmax2d.1.0e+03/mfpPerpCompos_slab%1.2f_sig%1.2f_Lc2D.%1.1e_LcSlab.%1.1e.dat' % (dir_theory, perc_slab, sig, Lc_2d, Lc_slab)
fname_theory_2	= '%s/xminSlab.1.0e-08_xmaxSlab1.0e+07_xmin2d.1.0e-06_xmax2d.1.0e+05/mfpPerpCompos_slab%1.2f_sig%1.2f_Lc2D.%1.1e_LcSlab.%1.1e.dat' % (dir_theory, perc_slab, sig, Lc_2d, Lc_slab)
fname_shalchi04_sim  = '../../../literatura/shalchi_etal04_fig5.png.dat'
fname_shalchi04_theo = '../../../literatura/shalchi_etal04_fig5.png_theory.dat'
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
theory_mfpPerp	= data_theory[1] # [AU]

#------------------------ theory WNLT (otro rango de integracion)
data_theory2	= loadtxt(fname_theory_2, unpack=True)
theory2_Rnorm	= data_theory2[0] # [1]
theory2_mfpPerp	= data_theory2[1] # [AU]

#------------------------ data del paper de Shalchi etal 2004
shal04_sim      = loadtxt(fname_shalchi04_sim, unpack=True)
shal04_theo     = loadtxt(fname_shalchi04_theo, unpack=True)


AUincm  = 1.5e13                        # [cm]
fig     = figure(1, figsize=(6, 4))
ax	    = fig.add_subplot(111)
axx	    = ax.twinx()

adim_RLs        = RLs/(AUincm*Lc_slab)
adim_MFPs_perp  = MFPs_perp/(AUincm*Lc_slab)

XLIMS	= array([3e-3, 2e0])
YLIMS	= array([1e-2, 1e0])

ax.plot(adim_RLs, adim_MFPs_perp, '-o', c='black', label='simulations')
ax.plot(theory_Rnorm, theory_mfpPerp/Lc_slab, '-', c='blue', label='WNLT approx', lw=3, alpha=.5)
ax.plot(theory2_Rnorm, theory2_mfpPerp/Lc_slab, 'v', c='blue', label='WNLT2 approx', lw=3, alpha=.5)
axx.plot(theory_Rnorm, theory_mfpPerp        , '-', c='red', lw=1, alpha=.0)
ax.plot(shal04_sim[0], shal04_sim[1]/Lc_slab, 'o', c='white', label='Sh, 2004')
ax.plot(shal04_theo[0], shal04_theo[1]/Lc_slab, '--', c='black')

TITLE1  = '$L_c^{2D}[AU]: %1.1e$   $L_c^{slab}[AU]: %1.1e$' % (Lc_2d, Lc_slab)
TITLE0  = '$Nm:%d$    $\sigma_{slab}^2/\sigma^2$: $%1.2f$    $(\sigma/Bo)^2$: $%1.1e$' %  (Nm, perc_slab, sig)
TITLE2 	= '$L_c^{2D} / L_c^{slab} = 0.1$'
TITLE   = TITLE0 + '\n' + TITLE2 + '\n'
ax.set_title(TITLE, fontsize=17.)
ax.set_xlabel('$R_L / L_c^{slab}$')
ax.set_ylabel('$\lambda_{\perp} / L_c^{slab}$')
axx.set_ylabel('[AU]')
ax.set_xscale('log'); axx.set_xscale('log')
ax.set_yscale('log'); axx.set_yscale('log')
ax.set_xlim(XLIMS);
ax.set_ylim(YLIMS); axx.set_ylim(YLIMS*Lc_slab)
ax.legend(loc='upper left')
ax.grid()

fname_fig       = '%s/mfpPerp.vs.Rnorm_Lc2d.%1.1e_LcSlab.%1.1e_sig.%1.2f_Nm%03d_slab.%1.2f_ii.png' % (dir_figs, Lc_2d, Lc_slab, sig, Nm, perc_slab)
savefig(fname_fig, format='png', dpi=135, bbox_inches='tight')
print " ----> generamos: %s" % fname_fig
close()
##
