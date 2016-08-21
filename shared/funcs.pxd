# distutils: language = c++
# Author: Jimmy J. Masias-Meza
import numpy as np
cimport numpy as np

# Numpy must be initialized. When using numpy from C or
# Cython you must _always_ do that, or you will have segfaults.
np.import_array()
#--- some datatypes 
#ctypedef np.int_t       Int
ctypedef np.float32_t   Float32
ctypedef np.ndarray     NDarray
#--- para q compile templates especificos
ctypedef double Doub
ctypedef int Int
