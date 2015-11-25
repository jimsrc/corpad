##
from funcs import *
import os

def generate_k_vs_t(Ek, dir_data):
	dir_info	= '%s/info' % dir_data
	fname_orient 	= '%s/orientations.in' % dir_info
	fname_plas	= '%s/plas.in' % dir_info
	fname_turb	= '%s/turb.in' % dir_info
	NPLAS		= count_lines_in_file(fname_orient)
	# orden del nro max de giroperiodos
	order_tmax	= int(log10(value(fname_plas, 'nmax_gyroperiods')))
	# nro de filas x archivo (nro de puntos q le pedi a la simulacion)
	nfil		= int(order_tmax*value(fname_plas, 'npoints') + 1)
	ncol		= 6			# nro de columnas x archivo: (t, x, y, z, mu, err-gamma)
	NBs		= int(value(fname_plas, 'nro_Bfield_realizations'))
	rigidity	= value(fname_plas, 'rigidity')
	#--------------
	Nm		= int(value(fname_turb, 'n_modos'))
	Bo		= value(fname_turb, 'Bunif')
	Sig		= value(fname_turb, 'sigma_Bo_ratio')
	perc_2d		= value(fname_turb, 'percent_2d')
	perc_slab	= value(fname_turb, 'percent_slab')
	Lc_2d		= value(fname_turb, 'Lc_2d')
	Lc_slab		= value(fname_turb, 'Lc_slab')
	lambda_min	= value(fname_turb, 'lambda_min')
	lambda_max	= value(fname_turb, 'lambda_max')

	print " ------> Ek [eV]: %g" % Ek
	calc_k_versus_t(dir_data, Ek, Sig, NPLAS, NBs, nfil, ncol, Bo, 
			Lc_2d, Lc_slab, Nm, perc_slab)

def calc_k_versus_t(dir_data, Ek, Sig, NPLAS, NBs, nfil, ncol, Bo,
		Lc_2d, Lc_slab, Nm, perc_slab):
	dir_plots = '../../plots'
	#dir_data= '../../output/Ek.%1.1eeV/sig%d' % (Ek, Sig)
	"""dir_out	  = '../../post/Ek.%1.1eeV/Nm%03d/slab%1.2f/sig.%1.1e/Lc2D.%1.1e_LcSlab.%1.1e' % (Ek, Nm, perc_slab, Sig, Lc_2d, Lc_slab)
	try: os.system('mkdir %s' % dir_out)
	except: print ""
	"""
	dir_out	  = '../../post'
	fname_out = '%s/k_vs_t_Ek.%1.1eeV_Nm%03d_slab%1.2f_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e.dat' % (dir_out, Ek, Nm, perc_slab, Sig, Lc_2d, Lc_slab)
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
	print ""
	savetxt(fname_out, data_out, fmt='%12.2f')
