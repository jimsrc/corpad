from funcs import *
"""
Puedo comparar estos resultdos con los de Qin. Ver figura 3.4 de Shalchi 2009 (libro).
Ahi muestra k(t). Son perfiles muy parecidos.
"""
Ek		= 1e6		# [eV]
Bo		= 5e-5
fname_inp	= '../post/k_vs_t_Ek.%1.2eeV.dat' % Ek
dir_plots	= '../plots'
t, kxx, kyy, kzz = loadtxt(fname_inp, unpack=True)
wc      = calc_omega(Bo, Ek)

#-------------------------------------------- fit k_perp
cc	= t>2e3
b	= 1.8e19/wc
m	= 3.6e22/wc
xo	= 0.
px	= make_fit_kxx([(t/wc)[cc], kxx[cc]], [b, m, xo])
py	= make_fit_kxx([(t/wc)[cc], kyy[cc]], [b, m, xo])
kxx_fit	= px[0]
kyy_fit	= py[0]
#-------------------------------------------- fit k_parall
ccz	= t>0.
cc	= t>2e3
b	= 3.5e22*(.2/wc) # 7e21 # 3.5e22 #1.8e22/wc
m	= -7.4e24*(.2/wc) # -2e24 # -7.4e24 #-3.6e21/wc
xo	= 0.
pz	= make_fit_kxx([(t/wc)[ccz], kzz[ccz]], [b, m, xo])
kzz_fit	= pz[0]
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
#--------------------------------------------
fig2	= figure(2, figsize=(6, 4))
ax2	= fig2.add_subplot(111)

ax2.scatter(t, kzz, edgecolor='none', c='black', alpha=.3)
ax2.plot(t[ccz], fun_hyperbola(pz[0], pz[1], pz[2], t[ccz]/wc), c='green')

TITLE	= 'fit: y=b+m/(x-xo)\nKzz: %1.1e' % (kzz_fit)
ax2.set_title(TITLE)
ax2.set_xlabel('time [$\Omega^{-1}$]')
ax2.set_ylabel('[cm2/s]')
ax2.grid()
#----------------------------------------

fname_fig_perp		= "%s/kperp_asymptotic.fit_Ek.%1.2eeV.png" % (dir_plots, Ek)
fname_fig_parall	= "%s/kparall_asymptotic.fit_Ek.%1.2eeV.png" % (dir_plots, Ek)
fig1.savefig(fname_fig_perp, format='png', dpi=135, bbox_inches='tight')
fig2.savefig(fname_fig_parall, format='png', dpi=135, bbox_inches='tight')
close(fig1); close(fig2)
#show(); close()
