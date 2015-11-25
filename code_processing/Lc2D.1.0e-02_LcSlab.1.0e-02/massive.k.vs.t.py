from calc_k_vs_t import *

Ek		= 1e9	# [eV]
Nm		= 128
perc_slab	= 0.20
sig		= 1e0
Lc_2d		= 1e-2
Lc_slab		= 1e-2
dir_data	= '../../output/Ek.%1.1eeV/Nm%03d/slab%1.2f/sig.%1.1e/Lc2D.%1.1e_LcSlab.%1.1e' % (Ek, Nm perc_slab, sig, Lc_2d, Lc_slab)

generate_k_vs_t(Ek, dir_data)
