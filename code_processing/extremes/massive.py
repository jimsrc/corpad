from calc_k_vs_t import *

#Ek	= 1e7			# [eV]
NPLAS   = 242
NBs	= 25
nfil	= 81			# nro de filas x archivo (nro de puntos q le pedi a la simulacion)
ncol	= 6			# nro de columnas x archivo: (t, x, y, z, mu, err-gamma)
Bo	= 5e-5			# [G] campo guia
Sigmas	= [1e+2]
Nms	= [512]
	
Ek = 1e6
print " ------> Ek [eV]: %g" % Ek
for Sig in Sigmas:
	for Nm in Nms:
		calc_k_versus_t(Ek, Sig, NPLAS, NBs, nfil, ncol, Bo, Nm)
