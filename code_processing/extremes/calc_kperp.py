from funcs import *
"""
Puedo comparar estos resultdos con los de Qin. Ver figura 3.4 de Shalchi 2009 (libro).
Ahi muestra k(t). Son perfiles muy parecidos.
"""
Ek		= 1e6		# [eV]
Bo		= 5e-5
Sig		= 1e2
Nm		= 512
fname_inp	= '../../post/sig.%.0e/Nm%d/k_vs_t_Ek.1e%deV.dat' % (Sig, Nm, log10(Ek))
dir_plots	= '../../plots'
t, kxx, kyy, kzz = loadtxt(fname_inp, unpack=True)
wc      = calc_omega(Bo, Ek)

#-------------------------------------------- fit k_perp
cc	= t>300
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
TITLE	= 'fit: y=b+m/(x-xo)\nKxx: %1.1e   Kyy:%1.1e\nKperp: %1.2e' % (kxx_fit, kyy_fit, kperp)
ax1.set_title(TITLE)
ax1.set_xlabel('time [$\Omega^{-1}$]')
ax1.set_ylabel('[cm2/s]')
ax1.grid()
ax1.legend(loc='upper right')
#--------------------------------------------
fname_fig_perp		= "%s/kperp_asymptotic.fit_Ek.%1.2eeV_sig%d_Nm%d.png" % (dir_plots, Ek, Sig, Nm)
print " ---> generando: %s" % fname_fig_perp
fig1.savefig(fname_fig_perp, format='png', dpi=135, bbox_inches='tight')
close(fig1)
#show(); close()
