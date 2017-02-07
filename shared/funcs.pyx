# Author: Jimmy J. Masias-Meza
from libc.stdlib cimport free, malloc, calloc
from cpython cimport PyObject, Py_INCREF#, PyMem_Malloc, PyMem_Free
from cpython.mem cimport PyMem_Malloc, PyMem_Free
#from cython.Utility.MemoryView import PyMem_New, PyMeM_Del # why doesn't work??
from cython.operator cimport dereference as deref
from libc.math cimport sqrt, sin, cos, M_PI, pow
from Bparker.Bparker import return_B as Bparker_vector
from numpy.linalg import norm
from numpy import array, zeros, float32 # no va en el .pxd!!

# agregamos la clase wrapper
#include "array_wrapper.pyx"

cpdef double calc_Rlarmor(Doub rigidity, Doub Bo):
    """
    input:
    Ek      : [eV] kinetic energy
    rigi..  : [V] rigidity
    Bo      : [G] magnetic field in Gauss
    output:
    Rl  : [cm] larmor radii
    """
    cdef:
        double q = (4.8032*1e-10) # [statC] carga PROTON
        double mo = 1.6726e-24 # [gr] masa PROTON
        double c = 3e10            # [cm/s] light speed
        double AU_in_cm = 1.5e13     # [cm]
        double E_reposo=938272013.0  # [eV] PROTON
        double beta, gamma, omg, v

    #rigidity = sqrt(Ek*Ek + 2.*Ek*E_reposo);
    #------------------------CALCULO DE GAMMA Y BETA
    gamma = pow(pow(rigidity/E_reposo,2) + 1. , 0.5)
    beta = pow(1. - 1/(gamma*gamma) , 0.5)
    #------------------------------CALCULO CICLOTRON
    omg = q * Bo / (gamma * mo * c)     # [s^-1]
    #---------------------------CALCULO RADIO LARMOR
    v   = beta * c              # [cm/s]
    #Rl[0]  = (v / omg) /AU_in_cm  # [AU]
    return (v / omg) # [cm]


def Bo_parker(r=1.0, th=M_PI/2., ph=0.0):
    """ input:
    r  [AU]  : heliodistance
    th [rad] : spherical polar angle (co-latitude)
    ph [rad] : azimuth
    default: (r=1.0, th=pi/2, ph=0.0)
    """
    #--- calc B-field at position (r, th)
    pos = array([r,th,ph], dtype=float32) #[AU], [rad], [rad]: position
    Bo  = norm(Bparker_vector(pos)) # [Gauss]
    return Bo

def Lc_memilia(r=1.0):
    """ 
    formula from Maria Emilia's thesis (Sec 4.4, p.88).
    """
    #Lc = 0.89*(r**(0.43))*1e6/(1e8) # [AU]
    Lc = 0.0059*(r**(0.43)) # [AU]
    return Lc

#EOF
