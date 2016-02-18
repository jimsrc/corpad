from funcs import *
"""
Puedo comparar estos resultdos con los de Qin. Ver figura 3.4 de Shalchi 2009 (libro).
Ahi muestra k(t). Son perfiles muy parecidos.
"""
def fit_kperp_asintot(Ek, dir_data, t_decr):
	dir_info        = '%s/info' % dir_data
	fname_plas      = '%s/plas.in' % dir_info
	fname_turb      = '%s/turb.in' % dir_info
	fname_orient    = '%s/orientations.in' % dir_info
	NPLAS           = count_lines_in_file(fname_orient)
	order_tmax      = int(log10(value(fname_plas, 'nmax_gyroperiods')))
	nfil            = int(order_tmax*value(fname_plas, 'npoints') + 1)
	ncol            = 6

	Nm              = int(value(fname_turb, 'n_modos'))
	Bo              = value(fname_turb, 'Bunif')
	Sig             = value(fname_turb, 'sigma_Bo_ratio')
	perc_2d         = value(fname_turb, 'percent_2d')
	perc_slab       = value(fname_turb, 'percent_slab')
	Lc_2d           = value(fname_turb, 'Lc_2d')
	Lc_slab         = value(fname_turb, 'Lc_slab')
	lambda_min      = value(fname_turb, 'lambda_min')
	lambda_max      = value(fname_turb, 'lambda_max')

	dir_out		= '../../post'
	fname_inp	= '%s/k_vs_t_Ek.%1.1eeV_Nm%03d_slab%1.2f_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e.dat' % (dir_out, Ek, Nm, perc_slab, Sig, Lc_2d, Lc_slab)
	dir_plots	= '../../plots'
	t, kxx, kyy, kzz = loadtxt(fname_inp, unpack=True)
	wc      = calc_omega(Bo, Ek)

	#-------------------------------------------- fit k_perp
	#cc	= t>5000
	cc	= t>t_decr
	b	= 1.8e19/wc
	m	= 3.6e22/wc
	xo	= 0.
	px	= make_fit_kxx([(t/wc)[cc], kxx[cc]], [b, m, xo])
	py	= make_fit_kxx([(t/wc)[cc], kyy[cc]], [b, m, xo])
	kxx_fit	= px[0]
	kyy_fit	= py[0]
	#--------------------------------------------
	fig1	= figure(1, figsize=(6, 4))
	ax1	= fig1.add_subplot(111)

	ax1.scatter(t, kxx, edgecolor='none', c='red', alpha=.3)
	ax1.scatter(t, kyy, edgecolor='none', c='blue', alpha=.3)
	ax1.scatter(t[cc], kxx[cc], edgecolor='none', c='red', label='kxx')
	ax1.scatter(t[cc], kyy[cc], edgecolor='none', c='blue', label='kyy')
	ax1.plot(t[cc], fun_hyperbola(px[0], px[1], px[2], t[cc]/wc), c='red')
	ax1.plot(t[cc], fun_hyperbola(py[0], py[1], py[2], t[cc]/wc), c='blue')

	kperp	= (kxx_fit + kyy_fit)/2.
	FIT_RESULTS = 'fit: y=b+m/(x-xo)\nKxx: %1.1e   Kyy:%1.1e\nKperp: %1.2e' % (kxx_fit, kyy_fit, kperp)
	SIMULAT_PARAMS	= 'Ek[eV]:%1.1e    perc_slab:%1.2f\nNm:%d    $(\sigma/Bo)^2$:%1.1e\n$Lc^{2D}$[AU]:%1.1e    $Lc^{slab}$[AU]:%1.1e' % (Ek, perc_slab, Nm, Sig, Lc_2d, Lc_slab)
	TITLE	= '%s\n%s' % (SIMULAT_PARAMS, FIT_RESULTS)
	ax1.set_title(TITLE)
	ax1.set_xlabel('time [$\Omega^{-1}$]')
	ax1.set_ylabel('[cm2/s]')
	ax1.grid()
	#--------------------------------------------
	fname_fig_perp		= "%s/kperp_asymptotic.fit_Ek.%1.2eeV_Nm%03d_slab%1.2f_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e.png" % (dir_plots, Ek, Nm, perc_slab, Sig, Lc_2d, Lc_slab)
	print " ---> generando: %s" % fname_fig_perp
	fig1.savefig(fname_fig_perp, format='png', dpi=135, bbox_inches='tight')
	close(fig1)
	#show(); close()
	return kperp

def fit_kparall_asintot(Ek, dir_data, t_decr):
	dir_info        = '%s/info' % dir_data
	fname_plas      = '%s/plas.in' % dir_info
	fname_turb      = '%s/turb.in' % dir_info
	fname_orient    = '%s/orientations.in' % dir_info
	NPLAS           = count_lines_in_file(fname_orient)
	order_tmax      = int(log10(value(fname_plas, 'nmax_gyroperiods')))
	nfil            = int(order_tmax*value(fname_plas, 'npoints') + 1)
	ncol            = 6

	Nm              = int(value(fname_turb, 'n_modos'))
	Bo              = value(fname_turb, 'Bunif')
	Sig             = value(fname_turb, 'sigma_Bo_ratio')
	perc_2d         = value(fname_turb, 'percent_2d')
	perc_slab       = value(fname_turb, 'percent_slab')
	Lc_2d           = value(fname_turb, 'Lc_2d')
	Lc_slab         = value(fname_turb, 'Lc_slab')
	lambda_min      = value(fname_turb, 'lambda_min')
	lambda_max      = value(fname_turb, 'lambda_max')

	dir_out		= '../../post'
	fname_inp	= '%s/k_vs_t_Ek.%1.1eeV_Nm%03d_slab%1.2f_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e.dat' % (dir_out, Ek, Nm, perc_slab, Sig, Lc_2d, Lc_slab)
	dir_plots	= '../../plots'
	t, kxx, kyy, kzz = loadtxt(fname_inp, unpack=True)
	wc      = calc_omega(Bo, Ek)

	#-------------------------------------------- fit k_perp
	#cc	= t>5000
	cc	= t>t_decr
	b	= 3.5e22*(.2/wc)
	m	= -7.4e24*(.2/wc)
	xo	= 0.
	pz      = make_fit_kxx([(t/wc)[cc], kzz[cc]], [b, m, xo])
	kzz_fit = pz[0]
	#--------------------------------------------
	fig1	= figure(1, figsize=(6, 4))
	ax1	= fig1.add_subplot(111)

	ax1.scatter(t, kzz, edgecolor='none', c='black', alpha=.3)
	ax1.scatter(t[cc], kzz[cc], edgecolor='none', c='black', label='kzz')
	ax1.plot(t[cc], fun_hyperbola(pz[0], pz[1], pz[2], t[cc]/wc), c='green', linewidth=3, alpha=.6)

	FIT_RESULTS = 'fit: y=b+m/(x-xo)    Kzz: %1.1e' % (kzz_fit)
	SIMULAT_PARAMS	= 'Ek[eV]:%1.1e    perc_slab:%1.2f\nNm:%d    $(\sigma/Bo)^2$:%1.1e\n$Lc^{2D}$[AU]:%1.1e    $Lc^{slab}$[AU]:%1.1e' % (Ek, perc_slab, Nm, Sig, Lc_2d, Lc_slab)
	TITLE	= '%s\n%s' % (SIMULAT_PARAMS, FIT_RESULTS)
	ax1.set_title(TITLE)
	ax1.set_xlabel('time [$\Omega^{-1}$]')
	ax1.set_ylabel('[cm2/s]')
	ax1.grid()
	#--------------------------------------------
	fname_fig_zz		= "%s/kzz_asymptotic.fit_Ek.%1.2eeV_Nm%03d_slab%1.2f_sig.%1.1e_Lc2d.%1.1e_LcSlab.%1.1e.png" % (dir_plots, Ek, Nm, perc_slab, Sig, Lc_2d, Lc_slab)
	print " ---> generando: %s" % fname_fig_zz
	fig1.savefig(fname_fig_zz, format='png', dpi=135, bbox_inches='tight')
	close(fig1)
	#show(); close()
	return kzz_fit
