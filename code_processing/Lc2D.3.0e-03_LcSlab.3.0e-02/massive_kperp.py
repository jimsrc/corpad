from calc_kperp import *
from funcs import *
import console_colors as ccl

Nm		= 128
sig		= 1e0
Lc_2d		= 3e-3
Lc_slab		= 3e-2

Ek	= 6e5				# [eV]
SLABs	= array([0.00, 0.01, 0.02, 0.05, 0.10, 0.20, 0.40, 0.60, 1.00])
n	= len(SLABs)
kperp	= zeros(n)
kzz	= zeros(n)

for i in range(n):
	perc_slab	= SLABs[i]
	dir_data        = '../../output/Ek.%1.1eeV/Nm%03d/slab%1.2f/sig.%1.1e/Lc2D.%1.1e_LcSlab.%1.1e' % (Ek, Nm, perc_slab, sig, Lc_2d, Lc_slab)
	t_decr		= 800 #TDECRs[i]  #200#200 #500 # 1000
	kperp[i] 	= fit_kperp_asintot(Ek, dir_data, t_decr)
	kzz[i] 		= fit_kparall_asintot(Ek, dir_data, 0.)

#---------- Mean Free Paths:
MFPs_perp	= calc_mfp(kperp, Ek)	# [cm]
MFPs_parall	= calc_mfp(kzz, Ek)	# [cm]

#---------- guardamos resultados en ascii
Bo		= 5e-5
RLs		= larmor(Bo, Ek)
data_out	= array([SLABs, kperp, kzz, MFPs_perp, MFPs_parall])
dir_out		= '../../post'
fname_out	= '%s/diff.pars_Nm%03d_Ek.%1.1e_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e.dat' % (dir_out, Nm, Ek, sig, Lc_2d, Lc_slab)
savetxt(fname_out, data_out.T, fmt='%5.5e')
print ccl.Rn + " ------> hemos creado: " + fname_out + ccl.W
