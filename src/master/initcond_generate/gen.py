from pylab import *
from numpy import *

Nth	= 16
Nphi	= 8
dth	= pi/Nth
dphi	= 2.*pi/Nphi

fname_out = "../../inputs/orientations_isotropic_Nth%d_Nph%d.in" % (Nth, Nphi)

data	= []
k	= 0
for i in range(Nth):
	for j in range(Nphi):
		cond = ((i==0) & (j==0)) | (i>0)
		if cond:
			print i," ", j, " cond:", cond
			#pause(1)
			data += [[i*dth, j*dphi]]
			k +=1
data += [[pi, 0.0]]	# agrego la direccion hacia abajo (mu=-1)
data = array(data)
savetxt(fname_out, data)
