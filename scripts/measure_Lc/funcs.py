#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
#from numpy.linalg import norm
from src_Bmodel import Bmodel
from shared.funcs import calc_Rlarmor, Bo_parker, Lc_memilia


class LcMeasure(object):
    def __init__(self, pd):
        self.m = Bmodel.Bmodel()
        self.m._build_pturb(pd=pd)

    def BiBj_profiles(self, ro, dr):
        """
        dada una realizacion de turbulencia, medir los 
        productos Bi(ro)*Bj(ro+dr) para diferentes valores `dr`.
        Output:
        - array of such product.
        """
        n = dr.size
        B1x, B2x = np.zeros((2,n), dtype=np.float32)
        B1y, B2y = np.zeros((2,n), dtype=np.float32)
      
        # `[:2]` is because I only want the fluctuating 
        # components (i.e. x and y components) of B.
        B1x[:], B1y[:] = self.m.Bxyz(ro)[:2] # evaluate at point "1"
        for i in range(n):
            B2x[i], B2y[i] = self.m.Bxyz(ro+dr[i])[:2] # eval at point "2"

        # B1 and B2 prefixes stand for "field evaluated and point 1
        # and point 2". So `B1x*B2x` is the correlation function
        # evaluated at the distance difference between point 1 and 
        # point 2. Actually `B1x*B2x` is an array as a function of an
        # array of difference-distances, using several points 2 and only
        # one point 1.
        _Rxx = np.prod([B1x,B2x], axis=0)
        _Ryy = np.prod([B1y,B2y], axis=0)
        return _Rxx, _Ryy

    def one_R_realiz(self, Nro=20, dr=None, nB=0, direcc='perp'):
        """
        Nro : nro of test-values for `ro`
        dr  : resolution of profile
        nB  : nro of turbulence-realization
        direcc: in which direction we build the correlation-function profile
        """
        self.m._build_par(nB=nB)

        Lc = self.m.read_param('Lc_slab')
        Rxx = np.zeros(dr.size, dtype=np.float32)
        Ryy = np.zeros(dr.size, dtype=np.float32)
        # set of "origin points"b
        xo_set = np.arange(0., Nro*Lc, 0.5*Lc) 

        if direcc=='perp':
            for xo, i in zip(xo_set, range(Nro)):
                ro      = [xo,0.,0.]
                _rxx, _ryy = self.BiBj_profiles(ro, dr)
                Rxx[:] += _rxx
                Ryy[:] += _ryy

            for xo, i in zip(xo_set, range(Nro)):
                ro      = [0.,xo,0.]
                _rxx, _ryy = self.BiBj_profiles(ro, dr)
                Rxx[:] += _rxx
                Ryy[:] += _ryy

            # i want the average (over all `ro`)
            Rxx /= 2.*Nro
            Ryy /= 2.*Nro

        elif direcc=='parall':
            for xo, i in zip(xo_set, range(Nro)):
                ro      = [0.,0.,xo]
                _rxx, _ryy = self.BiBj_profiles(ro, dr)
                Rxx[:] += _rxx
                Ryy[:] += _ryy

            # i want the average (over all `ro`)
            Rxx /= 1.*Nro
            Ryy /= 1.*Nro

        return Rxx, Ryy

#EOF
