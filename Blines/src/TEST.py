import numpy as np
import Bline as bl

fnameSW = '../../inputs/INPUT_TURB.SW.inp'
fnameSH = '../../inputs/INPUT_TURB.SH.inp'
fnameMC = '../../inputs/INPUT_TURB.MC.inp'
fout	= './tt.txt'
DX	= 0.1			# [AU]
ds	= 1e-4			# [AU]
every	= 100			# [1]
npoints = 600			# [1] nro de ptos por linea de campo
#pmod = mm.MODEL_TURB(fname)

pm = bl.Bline()
pm.set_Bmodel(DX, fnameSW, fnameSH, fnameMC)
pm.set_Blineparams(ds, every, fout, npoints)
pm.report()


pm.calc_Bline([-.5, 0., 0.])	# [AU]
d = np.array(pm.returnMat())

np.savetxt('ttt.txt', d)
"""
b = pm.Bfield()
bb = np.array(b)

print " bb: ", bb
"""
