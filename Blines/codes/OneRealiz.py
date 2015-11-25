##
from funcs import *
import Bline as bl
#------------------------------------------------- PARTICULA
#Rl		= 1.929785E-04  # 1e6eV
#Rl		= 1.980009E-03  # 1e8eV
Ek		= 1e8
npla		= 4
fnamePLA	= '../../output/single_orbits/Ek.%1.1eeV/Nm128/slab0.20/sig.1.0e+00/Lc2D.1.0e-02_LcSlab.1.0e-02/B00_pla%03d_traj.dat' % (Ek, npla)
datapla		= loadtxt(fnamePLA, unpack=True)
Xpla		= datapla[1]
Ypla		= datapla[2]
Zpla		= datapla[3]


# inicia params del campo
fnameTURB 	= '../../inputs/single_orbits/TURB.sig1.0e+00_slab0.2_Bo.05.in'
fout		= './tt.txt'
slab_perc	= float(leer_param(fnameTURB, 'percent_slab'))
pm 		= bl.Bline()
pm.set_Bmodel(fnameTURB)

# inicia params de linea de campo
ds	= 0.1*pm.retLambdaMin()/AUincm 	#.5e-5			# [AU]
every	= 200				# [1]
npoints = 100				# [1] nro de ptos por linea de campo (*)
pm.set_Blineparams(ds, every, fout, npoints)
pm.report()
Bline_size = every*npoints*ds
sigma	= pm.retSigmaBoRatio()
# NOTA: de manera que npoints*every = nro total de ptos computados
#----------------------------------------------------
Lc	= 0.01		# [AU]
ny= 10; nx= 12
Xs	= linspace(-Lc, Lc, nx)
Ys	= linspace(-0.5*Lc, 0.5*Lc, ny)
zo	= -0.5*Lc

n_Brealiz = 1
n_Blines  = 10
TITLE0	  = "Bline_size[$L_c^{slab}$]: %g\n n Blines(nx, ny): %dx%d" % (Bline_size/Lc, nx, ny)

for k in range(n_Brealiz):
	print " ======== >>> B-realization %d" % k
	fig1 	= figure(1, figsize=(9, 7))
	ax	= fig1.add_subplot(1,1,1)
	ax.set_title(TITLE0)
	#for j in range(n_Blines):
	for yo in Ys:
		for xo in Xs:
			print " xo, yo: %f, %f" % (xo, yo)
			pm.calc_Bline([xo, yo, zo])
			d = array(pm.returnMat())
			x 	= d[:, 0]
			y 	= d[:, 1]
			z 	= d[:, 2]
			bmod	= d[:, 6]
			#------------------------------------
			ALPHA = 0.4*(yo-Ys[0])/(Ys[ny-1]-Ys[0]) + 0.1
			ax.plot(x/Lc, z/Lc, 'o-', c='black', markeredgecolor='none', alpha=ALPHA, markersize=1.)
			ax.plot(Xpla/Lc, Zpla/Lc, '-', c='blue', markeredgecolor='none', markersize=.7)
			ax.scatter(x[0]/Lc, z[0]/Lc, c='red', edgecolor='none', s=3.)
	#---------------------------------------
	ax.set_xlabel('$X / L_c^{slab}$')
	ax.set_ylabel('$Z / L_c^{slab}$')
	ax.grid()
	#---------------------------------------
	fname_out = 'prof_sig%1.2f_slab%1.2f_B%02d' % (sigma, slab_perc, k)
	fname_fig = '../plots/B00.realiz/%s_pla%03d.png' % (fname_out, npla)
	fname_odat= '../data/%s.txt' % fname_out
	#show()
	fig1.savefig(fname_fig, format='png', dpi=100, bbox_inches='tight')
	print " ----> se genero: %s" % fname_fig
	#close(fig0)
	close(fig1)
	savetxt(fname_odat, d)
	print " ----> se genero: %s" % fname_odat

	pm.next_Brealization()

