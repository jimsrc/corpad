from calc_k_vs_t import *
import os

Ek		= 6e5	# [eV]
Nm		= 128
perc_slab	= 0.01 #0.02 #0.05 #0.10 # 0.00, 0.20, 0.40, 0.60, 1.00
sig		= 1e0
Lc_2d		= 3e-3
Lc_slab		= 3e-2

Bo  = 5e-5
dir_data	= '../../../output/Ek.%1.1eeV/Nm%03d/slab%1.2f/sig.%1.1e/Lc2D.%1.1e_LcSlab.%1.1e' % (Ek, Nm, perc_slab, sig, Lc_2d, Lc_slab)

kt = k_vs_t(Ek, dir_data)
kt.calc_k_versus_t(Bo, './test')

#generate_k_vs_t(Ek, dir_data)
