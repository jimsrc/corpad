from calc_k_vs_t import *

#Ek	= 1e7			# [eV]
NPLAS   = 242
NBs	= 50
nfil	= 81			# nro de filas x archivo (nro de puntos q le pedi a la simulacion)
ncol	= 6			# nro de columnas x archivo: (t, x, y, z, mu, err-gamma)
Bo	= 5e-5			# [G] campo guia

Nm	= 4
for i in range(6, 7):
	Ek = pow(10, i)
	print " ------> Ek [eV]: %g" % Ek
	calc_k_versus_t(Ek, Nm, NPLAS, NBs, nfil, ncol, Bo)
