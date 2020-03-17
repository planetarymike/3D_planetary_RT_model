# py_corona_sim.pyx -- Cython wrapper for H corona simulation object
# distutils: language = c++
# cython: language_level=3


from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

# Import the Python-level symbols of numpy
import numpy as np
cimport numpy as cnp

# Import the C-level symbols of numpy
#cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
#np.import_array()

cdef extern from "observation_fit.h":
    cdef cppclass observation_fit:
        observation_fit() except +
        void add_observation(vector[vector[double]], vector[vector[double]])
        void set_g_factor(double)
        void generate_source_function(double, double)
        vector[double] brightness()
        
cdef class Pyobservation_fit:
    cdef observation_fit *thisptr #holds the reference to the cpp class
    def __cinit__(self):
        self.thisptr = new observation_fit()
    def __dealloc__(self):
        del self.thisptr
    def add_observation(self, double[:,:] loc_arr, double[:,:] dir_arr):
        cdef vector[vector[double]] loc_vec,dir_vec
        loc_vec.resize(loc_arr.shape[0])
        dir_vec.resize(dir_arr.shape[0])
        for i in range(loc_arr.shape[0]):
            loc_vec[i].resize(3)
            dir_vec[i].resize(3)
            for j in range(3):
                loc_vec[i][j] = loc_arr[i,j]
                dir_vec[i][j] = dir_arr[i,j]
        self.thisptr.add_observation(loc_vec,dir_vec)
    def set_g_factor(self, g):
        self.thisptr.set_g_factor(g)
    def generate_source_function(self, double nH, double T):
        self.thisptr.generate_source_function(nH,T)
    def brightness(self):
        return np.asarray(self.thisptr.brightness())
