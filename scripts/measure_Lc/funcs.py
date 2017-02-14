#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
from numpy import cos, sin
#from numpy.linalg import norm
from src_Bmodel import Bmodel
from shared.funcs import calc_Rlarmor, Bo_parker, Lc_memilia


class LcMeasure(object):
    def __init__(self, pd, bidim=False):
        self.m = Bmodel.Bmodel()
        self.m._build_pturb(pd=pd)
        self.bidim = bidim # True: make bidim plot array of R-funcs

    def BiBj_profiles(self, ro, dr):
        """
        dada una realizacion de turbulencia, medir los 
        productos Bi(ro)*Bj(ro+dr) para diferentes valores `dr`.
        Output:
        - array of such product.
        """
        n = dr[:,0].size # size along one of its coordinates
        B1x, B2x = np.zeros((2,n), dtype=np.float32)
        B1y, B2y = np.zeros((2,n), dtype=np.float32)
      
        # `[:2]` is because I only want the fluctuating 
        # components (i.e. x and y components) of B.
        B1x[:], B1y[:] = self.m.Bxyz(ro)[:2] # evaluate at point "1"
        for i in range(n):
            r = ro + dr[i,:] # this is point "2"
            B2x[i], B2y[i] = self.m.Bxyz(r)[:2] # eval at point "2"

        # B1 and B2 prefixes stand for "field evaluated and point 1
        # and point 2". So `B1x*B2x` is the correlation function
        # evaluated at the distance difference between point 1 and 
        # point 2. Actually `B1x*B2x` is an array as a function of an
        # array of difference-distances, using several points 2 and only
        # one point 1.
        _Rxx = np.prod([B1x,B2x], axis=0)
        _Ryy = np.prod([B1y,B2y], axis=0)
        return _Rxx, _Ryy

    def one_R_realiz(self, Nro=20, dr=None, dth=None, nB=0, direcc='perp'):
        """
        Nro : nro of test-values for `ro`
        dr  : r-resolution of profile (r is the radii)
        nB  : nro of turbulence-realization
        direcc: in which direction we build the correlation-function profile
        """
        self.m._build_par(nB=nB)

        Lc = self.m.read_param('Lc_slab')

        # set of "origin points"b
        # NOTE: we move the origin in a scale comparable to `Lc`.
        xo_set = np.arange(0., Nro*Lc, 0.5*Lc) 

        #--- bidimensional calculus for R-funcs
        if self.bidim:
            Rxx = np.zeros((dr.size,dth.size), dtype=np.float32)
            Ryy = np.zeros((dr.size,dth.size), dtype=np.float32)

            _dph = 0.0 # azimuth (in plane x-y) # TODO: loop 
            _dr = np.empty((dr.size,3), dtype=np.float32)

            # loop over the origin-positions
            for xo, i in zip(xo_set, range(Nro)):
                for ith in range(dth.size): # loop over theta (polar angle)
                    _dth = dth[ith] # polar angle (angle respect to z-axis)
                    _dr[:,0] = dr[:]*cos(_dph)*sin(_dth)
                    _dr[:,1] = dr[:]*sin(_dph)*sin(_dth)
                    _dr[:,2] = dr[:]*cos(_dth)
                    ro       = [xo,0.,0.] # TODO: randomizar direccion del origen 'ro'
                    _rxx, _ryy = self.BiBj_profiles(ro, _dr)
                    Rxx[:,ith] += _rxx
                    Ryy[:,ith] += _ryy

            Rxx /= 1.*xo_set.size
            Ryy /= 1.*xo_set.size

            return Rxx, Ryy

        Rxx = np.zeros(dr.size, dtype=np.float32)
        Ryy = np.zeros(dr.size, dtype=np.float32)
        # `dr` is assumed in the x-y plane
        if direcc=='perp':
            _dr = np.empty((dr.size,3), dtype=np.float32)

            #--- along x
            _dr[:,0] = dr # displacements are on x-direction
            _dr[:,1] = 0.
            _dr[:,2] = 0.
            for xo, i in zip(xo_set, range(Nro)):
                ro      = [xo,0.,0.]
                _rxx, _ryy = self.BiBj_profiles(ro, _dr)
                Rxx[:] += _rxx
                Ryy[:] += _ryy

            #--- along y
            _dr[:,0] = 0.
            _dr[:,1] = dr # displacements are on y-direction
            _dr[:,2] = 0.
            for xo, i in zip(xo_set, range(Nro)):
                ro      = [0.,xo,0.]
                _rxx, _ryy = self.BiBj_profiles(ro, _dr)
                Rxx[:] += _rxx
                Ryy[:] += _ryy

            # i want the average (over all `ro`)
            Rxx /= 2.*Nro
            Ryy /= 2.*Nro

        # `dr` is assumed in the z-direction
        elif direcc=='parall':
            _dr = np.empty((dr.size,3), dtype=np.float32)
            _dr[:,0] = 0.
            _dr[:,1] = 0.
            _dr[:,2] = dr
            for xo, i in zip(xo_set, range(Nro)):
                ro      = [0.,0.,xo]
                _rxx, _ryy = self.BiBj_profiles(ro, _dr)
                Rxx[:] += _rxx
                Ryy[:] += _ryy

            # i want the average (over all `ro`)
            Rxx /= 1.*Nro
            Ryy /= 1.*Nro

        else:
            raise SystemError("""
            direcc argument must be one of these: 
            - perp
            - parall
            - axisymmet
            """)

        return Rxx, Ryy

#EOF
