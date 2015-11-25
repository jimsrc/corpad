from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
AUincm	= 1.5e13		# [cm]

def leer_param(fname, varname):
	f = open(fname, 'r')
	for line in f:
		ll = line.split()
		if(ll[1]==varname):
			return str(ll[0])
	print " varname '%s' NOT FOUND in; %s" % (varname, fname)
	f.close()
	return 'nan'
