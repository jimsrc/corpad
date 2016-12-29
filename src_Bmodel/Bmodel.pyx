# distutils: language = c++
# Author: Jimmy J.
#from libcpp.string cimport string

# Declare the prototype of the C function we are interested in calling

from libc.stdlib cimport free, malloc, calloc
from cpython cimport PyObject, Py_INCREF#, PyMem_Malloc, PyMem_Free
from cpython.mem cimport PyMem_Malloc, PyMem_Free
#from cython.Utility.MemoryView import PyMem_New, PyMeM_Del # why doesn't work??
from cython.operator cimport dereference as deref
from libc.math cimport sqrt, sin, cos

# agregamos la clase wrapper
include "array_wrapper.pyx"


cdef class Bmodel(object):
    cdef PARAMS_TURB            *pt
    cdef MODEL_TURB             *mt

    def __cinit__(self):
        pass

    def _build_pturb(s, dict pd=None):
        """ build PARAMS_TURB object 'self.pt' from
        dictionary input 'self.pdict' """
        s.pt = new PARAMS_TURB()
        if s.pt is NULL:
            raise MemoryError()

        #pd = s.pdict
        # parametros fisicos
        s.pt.Nm_slab = pd['Nm_slab']
        s.pt.Nm_2d   = pd['Nm_2d']
        s.pt.lmin_s  = pd['lmin_s']
        s.pt.lmax_s  = pd['lmax_s']
        s.pt.lmin_2d = pd['lmin_2d']
        s.pt.lmax_2d = pd['lmax_2d']
        s.pt.Lc_slab = pd['Lc_slab']
        s.pt.Lc_2d   = pd['xi']*pd['Lc_slab']
        s.pt.sigma_Bo_ratio = pd['sigma_Bo_ratio']
        s.pt.percent_slab = pd['ratio_slab']
        s.pt.percent_2d   = 1.0-pd['ratio_slab'] #pd['percent_2d']
        #s.pt.Bo      = pd['Bo']
        # semillas
        s.pt.sem.slab[0] = pd['sem_slab0']
        s.pt.sem.slab[1] = pd['sem_slab1']
        s.pt.sem.slab[2] = pd['sem_slab2']
        s.pt.sem.two[0]  = pd['sem_two0']
        s.pt.sem.two[1]  = pd['sem_two1']

    def _build_par(s, int nB=0):
        """ objeto MODEL_TURB (todo):
        - aloco memoria para dB, B, etc..
        - defino modelo B con 'self.pt'
        - build dB spectra
        - fix B realization
        """
        s.mt        = new MODEL_TURB()
        if s.mt is NULL:
            raise MemoryError()

        ndim        = 3
        #s.mt.B       = PyMem_New(double, 3)
        # use 'PyMem_Malloc' instead of 'malloc'
        # src: https://docs.python.org/2/c-api/memory.html
        s.mt.B       = <Doub*> calloc(ndim, sizeof(Doub))
        s.mt.dB      = <Doub*> calloc(ndim, sizeof(Doub))
        s.mt.dB_SLAB = <Doub*> calloc(ndim, sizeof(Doub))
        s.mt.dB_2D   = <Doub*> calloc(ndim, sizeof(Doub))

        #NOTE: al 's.pt' lo construi en 'self._build_pturb()'
        s.mt.p_turb  = s.pt[0] #le paso todos los parametros! (MAGIA!!??!?)
        s.mt.p_turb.build_spectra()
        s.mt.fix_B_realization(nB=nB)

    def read_param(self, name):
        """
        read parameters from the built B-turbulence-model (see
        _build_pturb && _build_par methods).
        """
        cdef double *ptr
        cdef np.ndarray ndarray

        if name=='Bk_SLAB':
            n = self.mt.p_turb.Nm_slab  # size of array
            ptr = &(self.mt.p_turb.Bk_SLAB[0]) # set the pointer
        elif name=='Bk_2D':
            n = self.mt.p_turb.Nm_2d
            ptr = &(self.mt.p_turb.Bk_2D[0])
        elif name=='k_s':
            n = self.mt.p_turb.Nm_slab  # size of array
            ptr = &(self.mt.p_turb.k_s[0]) # set the pointer
        elif name=='k_2d':
            n = self.mt.p_turb.Nm_2d  # size of array
            ptr = &(self.mt.p_turb.k_2d[0]) # set the pointer
        elif name=='dk_s':
            n = self.mt.p_turb.Nm_slab  # size of array
            ptr = &(self.mt.p_turb.dk_s[0]) # set the pointer
        elif name=='dk_2d':
            n = self.mt.p_turb.Nm_2d  # size of array
            ptr = &(self.mt.p_turb.dk_2d[0]) # set the pointer
        elif name=='Lc_slab':
            return self.mt.p_turb.Lc_slab
        elif name=='Lc_2d':
            return self.mt.p_turb.Lc_2d
        else:
            return None
    
        #--- wrap the C++ array
        arrw = ArrayWrapper()   # numpy-array wrapper
        arrw.set_data(n, <void*> ptr, survive=True)
        ndarray = np.array(arrw, copy=False)
        ndarray.base = <PyObject*> arrw
        Py_INCREF(arrw)
        return ndarray


    def __dealloc__(self):
        free(self.mt.B)
        free(self.mt.dB)
        free(self.mt.dB_SLAB)
        free(self.mt.dB_2D)
        del self.pt

    def Bxyz(self, xyz):
        cdef double pos[3]
        pos[0]=xyz[0]; pos[1]=xyz[1]; pos[2]=xyz[2]
        self.mt.calc_B(&(pos[0]))
        B = np.zeros(3)
        B[0]=self.mt.B[0]; B[1]=self.mt.B[1]; B[2]=self.mt.B[2]
        return B

#EOF
