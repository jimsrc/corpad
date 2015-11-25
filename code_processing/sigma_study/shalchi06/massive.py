from calc_k_vs_t import *

#Ek	= 1e7			# [eV]
dir_src = '../../../output'
dir_out	= '../../../post'
dir_plt = '../../../plots'

NPLAS   = 242
NBs	= 50
nfil	= 81			# nro de filas x archivo (nro de puntos q le pedi a la simulacion)
ncol	= 6			# nro de columnas x archivo: (t, x, y, z, mu, err-gamma)
Bo	= 5e-5			# [G] campo guia
Sigmas	= [1e+0]
Nms	= [128]
	
Ek = 6e7
slab_perc = [0.01]
Sig = 1.0; Nm=128
print " ------> Ek [eV]: %g" % Ek
for slab_p in slab_perc:
	print " ----> s:", slab_p
	calc_k_versus_t(dir_src, dir_out, dir_plt, 
				Ek, Sig, NPLAS, NBs, nfil, ncol, Bo, Nm, slab_p)
