##
from funcs import *
"""
Ek	= 1e7			# [eV]
NPLAS   = 211
NBs	= 50
nfil	= 2001			# nro de filas x archivo (nro de puntos q le pedi a la simulacion)
ncol	= 6			# nro de columnas x archivo
Bo	= 5e-5			# [G] campo guia
"""
def calc_k_versus_t(Ek, Sig, NPLAS, NBs, nfil, ncol, Bo):
	dir_plots = '../../plots'
	dir_data= '../../output/output_Ek.%3.1fMeV/sig%d' % (Ek, Sig)
	fname_out = '../../post/sig%d/k_vs_t_Ek.%3.1fMeV.dat' % (Sig, Ek)

	RANGE_X	= [-1e2, 4.1e4]
	RANGE	= [-.6, .6]
	#---------------------
	# nok	: nro de files q existe Y tienen data
	# nbad	: nro de files q solicite y no existen
	# time	: grilla temporal
	DATA, time, nok, nbad = load_trajectories(NBs, NPLAS, nfil, ncol, dir_data)

	print " nro de plas: ", NPLAS
	print " nro de B-realizations: ", NBs
	print " nro de ptos por trayectoria: %d\n" % nfil
	print " nro de archivos q existe c/data: %d/%d " % (nok, nok+nbad)
	print " nro de archivos q pedi y NO existen: %d/%d " % (nbad, nok+nbad)
	#---------------------
	fig	= figure(1, figsize=(6,4))
	ax	= fig.add_subplot(111)
	axx	= ax.twinx()

	every	= 1			# no en c/tiempo, sino cada 'every'
	tt, x2, y2, z2 = sqr_deviations(DATA, time, every)

	AUinkm	= 1.5e8
	AUincm	= AUinkm*1e5		# [cm]
	r2	= x2 + y2
	r2	= r2*AUincm**2		# [cm^2]
	x2	= x2*AUincm**2		# [cm^2]
	y2	= y2*AUincm**2		# [cm^2]
	z2	= z2*AUincm**2		# [cm^2]
	wc	= calc_omega(Bo, Ek) #4.781066E-01 #4.325188E-01 #Ek=1e8eV  #4.735689E-01 # Ek=1e7eV #4.781066E-01 # Ek=1e6eV
	print " wc[s-1]: ", wc
	tt	= tt*wc				# [1]
	#-------------------
	kxx	= x2/(2.*tt/wc)		# [cm2/s]
	kyy	= y2/(2.*tt/wc)		# [cm2/s]
	kzz	= z2/(2.*tt/wc)		# [cm2/s]
	#-- guarda data kxx(t)
	data_out	= array([tt, kxx, kyy, kzz]).T
	data_out	= data_out[1:]				# el 1er tiempo no lo guardo xq es division por zero 1/2t
	print " ---> guardando: %s" % fname_out
	savetxt(fname_out, data_out, fmt='%12.2f')
