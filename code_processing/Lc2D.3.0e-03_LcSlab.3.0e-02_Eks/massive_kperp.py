from calc_kperp import *
from funcs import *
import console_colors as ccl

Nm		= 128
sig		= 1e0
Lc_2d		= 3e-3
Lc_slab		= 3e-2

EKs		= array([6e5, 6e6, 1e8, 2e8])				# [eV]
perc_slab	= 0.20
n	= len(EKs)
kperp	= zeros(n)
kzz	= zeros(n)

for i in range(n):
	Ek		= EKs[i]
	dir_data        = '../../output/Ek.%1.1eeV/Nm%03d/slab%1.2f/sig.%1.1e/Lc2D.%1.1e_LcSlab.%1.1e' % (Ek, Nm, perc_slab, sig, Lc_2d, Lc_slab)
	t_decr		= 400 #800 #TDECRs[i]  #200#200 #500 # 1000
	kperp[i] 	= fit_kperp_asintot(Ek, dir_data, t_decr)
	kzz[i] 		= fit_kparall_asintot(Ek, dir_data, 0.)

#---------- Mean Free Paths:
MFPs_perp	= calc_mfp(kperp, EKs)	# [cm]
MFPs_parall	= calc_mfp(kzz, EKs)	# [cm]

#---------- guardamos resultados en ascii
Bo		= 5e-5
RLs		= larmor(Bo, EKs)
data_out	= array([EKs, RLs, kperp, kzz, MFPs_perp, MFPs_parall])
dir_out		= '../../post'
fname_out	= '%s/diff.pars_Nm%03d_slab.%1.2f_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e.dat' % (dir_out, Nm, perc_slab, sig, Lc_2d, Lc_slab)
savetxt(fname_out, data_out.T, fmt='%5.5e')
print ccl.Rn + " ------> hemos creado: " + fname_out + ccl.W
