import numpy as np
import read_Bmodel as mm

fname = '../../inputs/INPUT_TURB_01.inp'
#pmod = mm.MODEL_TURB(fname)

pm = mm.Bmodel()
pm.buildB(1e-3, fname, fname, fname)

pm.calc_B([44.,33.,1e10])

b = pm.Bfield()
bb = np.array(b)

print " bb: ", bb
