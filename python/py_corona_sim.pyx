# py_corona_sim.pyx -- Cython wrapper for H corona simulation object
# distutils: language = c++

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

# Import the Python-level symbols of numpy
import numpy as np

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
    def add_observation(self, loc_arr, dir_arr):
        self.thisptr.add_observation(loc_arr,dir_arr)
    def set_g_factor(self, g):
        self.thisptr.set_g_factor(g)
    def generate_source_function(self, double nH, double T):
        self.generate_source_function(nH,T)
    def brightness(self):
        return np.asarray(self.thisptr.brightness())
