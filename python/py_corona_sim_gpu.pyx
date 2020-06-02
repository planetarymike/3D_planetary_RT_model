# py_corona_sim.pyx -- Cython wrapper for H corona simulation object
# distutils: language = c++
# cython: language_level=3
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

# Import the Python-level symbols of numpy
import numpy as np
cimport numpy as np

ctypedef float Real #GPU (for high throughput)

cdef extern from "observation_fit.hpp":
    cpdef cppclass observation_fit:
        observation_fit()
        void add_observation(vector[vector[Real]] MSO_locations, vector[vector[Real]] MSO_directions)
        void set_g_factor(Real &g)
        void generate_source_function(Real nHexo, Real Texo)
        vector[Real] brightness()
        
cdef class Pyobservation_fit:
    cpdef observation_fit *thisptr #holds the reference to the cpp class
    def __cinit__(self):
        self.thisptr = new observation_fit()
    def __dealloc__(self):
        del self.thisptr
    def add_observation(self, loc_arr, dir_arr):
        cdef vector[vector[Real]] loc_vec,dir_vec
        loc_vec.resize(loc_arr.shape[0])
        dir_vec.resize(dir_arr.shape[0])
        for i in range(loc_arr.shape[0]):
            loc_vec[i].resize(3)
            dir_vec[i].resize(3)
            for j in range(3):
                loc_vec[i][j] = np.float32(loc_arr[i,j])
                dir_vec[i][j] = np.float32(dir_arr[i,j])
        self.thisptr.add_observation(loc_vec,dir_vec)
    def set_g_factor(self, Real g):
        self.thisptr.set_g_factor(g)
    def generate_source_function(self, Real nH, Real T):
        self.thisptr.generate_source_function(nH,T)
    def brightness(self):
        return np.asarray(self.thisptr.brightness())
