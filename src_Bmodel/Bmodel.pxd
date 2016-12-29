#--- librerias de c
#from libc.math cimport sqrt, sin, cos
from libc.math cimport sqrt, pow


#--- estructuras c++
cdef extern from "nr3.h":
    #cdef typename Doub # declarar lo q ya esta definido en general.h?
    #------ NRvector ------#
    cdef cppclass NRvector[T]:
        NRvector(int n) except +
        NRvector(int n, const T *array) except +
        NRvector(int n, const T &value) except +
        #NRvector & operator=(NRvector &rhs)
        T& operator[](int i)
        int size() const

#--- para q compile templates especificos
ctypedef double Doub
ctypedef int Int
ctypedef NRvector[Doub] VecDoub
ctypedef NRvector[Doub] VecDoub_IO
ctypedef NRvector[Doub] VecDoub_O
ctypedef const NRvector[Doub] VecDoub_I


#--- clases PARAMS_... y MODEL_TURB
cdef extern from "defs_turb.h":
    cpdef cppclass PARAMS_SEM:
        PARAMS_SEM()
        long slab[3]
        long two[2]
    
    cpdef cppclass PARAMS_TURB:
        PARAMS_TURB()
        void build_spectra()
        Int Nm_slab, Nm_2d
        Doub lmin_s, lmax_s
        Doub lmin_2d, lmax_2d
        #Int n_modos
        #Doub lambda_min
        #Doub lambda_max
        Doub Lc_slab, Lc_2d
        Doub sigma_Bo_ratio;
        Doub percent_slab, percent_2d;
        Doub gS, g2D # potencia espectral slab/2D
        Doub Bo    # campo uniforme
        Doub sigma_S    # intensidad slab
        Doub sigma_2D   # intensidad 2D
        PARAMS_SEM sem  # semillas
        Doub *Bk_SLAB, *Bk_2D   # fourier intensities
        Doub *k_s, *k_2d        # fourier modes
        Doub *dk_s, *dk_2d;
        Doub gS, g2D            # spectral indexes

    cpdef cppclass MODEL_TURB:
        MODEL_TURB()
        double *B
        double *dB
        double *dB_2D
        double *dB_SLAB
        PARAMS_TURB p_turb
        void calc_B(const double *)
        void fix_B_realization(const int nB)




cdef double AU_in_cm #= 1.5e13
#AU_in_cm = 1.5e13 # corre ok, pero no funciona

#EOF
