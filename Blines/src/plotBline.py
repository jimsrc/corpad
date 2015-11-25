##
from funcs import *
import Bline as bl
# inicia params del campo
fnameTURB = '../../inputs/INPUT_TURB.ALONE.inp'
fout	= './tt.txt'
pm = bl.Bline()
pm.set_Bmodel(fnameTURB)

# inicia params de linea de campo
ds	= 0.1*pm.retLambdaMin()/AUincm 	#.5e-4			# [AU]
every	= 200				# [1]
npoints = 200				# [1] nro de ptos por linea de campo (*)
pm.set_Blineparams(ds, every, fout, npoints)
pm.report()
Bline_size = every*npoints*ds
sigma	= pm.retSigmaBoRatio()
# NOTA: de manera que npoints*every = nro total de ptos computados

#----------------------------------------------------
#fig0 	= figure(1, figsize=(8, 6))
#ax0	= fig0.gca(projection='3d')

xo	= -.3
yo	= -1.		# [AU]
dy	= .01		# [AU]
zo	= 0.0

n_Brealiz = 10
n_Blines  = 10
TITLE0	  = "Bline_size[AU]: %g\n n Blines: %d" % (Bline_size, n_Blines)

for k in range(n_Brealiz):
	print " ======== >>> B-realization %d" % k
	fig1 	= figure(1, figsize=(12, 9))
	ax	= []
	ax	+= [fig1.add_subplot(2,2,1)]
	ax	+= [fig1.add_subplot(2,2,2)]
	ax	+= [fig1.add_subplot(2,2,3)]
	ax	+= [fig1.add_subplot(2,2,4)]
	ax[0].set_title(TITLE0)
	for j in range(n_Blines):
		y = yo + j*dy
		pm.calc_Bline([xo, y, zo])
		d = array(pm.returnMat())
		x 	= d[:, 0]
		y 	= d[:, 1]
		z 	= d[:, 2]
		bmod	= d[:, 6]
		#------------------------------------
		"""ax0.plot(x, y, z, alpha=.4, c='blue', lw=2)
		ax0.scatter(x[0], y[0], z[0], c='red')
		ax0.set_xlabel('X [AU]')
		ax0.set_ylabel('Y [AU]')
		ax0.set_zlabel('Z [AU]')"""
		#------------------------------------
		ax[0].plot(x, y, 'o-', c='blue', markeredgecolor='none', alpha=.4, markersize=1.)
		ax[0].scatter(x[0], y[0], c='red')
		#
		ax[1].plot(x, z, 'o-', c='blue', markeredgecolor='none', alpha=.4, markersize=1.)
		ax[1].scatter(x[0], z[0], c='red')
		#
		ax[2].plot(y, z, 'o-', c='blue', markeredgecolor='none', alpha=.4, markersize=1.)
		ax[2].scatter(y[0], z[0], c='red')
		#
		ax[3].plot(x, bmod/1e-5, 'o-', c='black', markeredgecolor='none', alpha=.4, markersize=2.)

	#---------------------------------------
	ax[0].set_xlabel('X [AU]')
	ax[0].set_ylabel('Y [AU]')
	ax[1].set_xlabel('X [AU]')
	ax[1].set_ylabel('Z [AU]')
	ax[2].set_xlabel('Y [AU]')
	ax[2].set_ylabel('Z [AU]')
	ax[3].set_xlabel('X [AU]')
	ax[3].set_ylabel('B [nT]')
	
	ax[3].set_ylim(0, 25)
	for l in range(4):
		ax[l].grid()
	#---------------------------------------
	fname_out = 'prof_sig%1.2f_xo=%1.2fAU_B%02d' % (sigma, xo, k)
	fname_fig = '../plots/%s.png' % fname_out
	fname_odat= '../data/%s.txt' % fname_out
	#show()
	fig1.savefig(fname_fig, format='png', dpi=300, bbox_inches='tight')
	print " ----> se genero: %s" % fname_fig
	#close(fig0)
	close(fig1)
	savetxt(fname_odat, d)
	print " ----> se genero: %s" % fname_odat

	pm.next_Brealization()

