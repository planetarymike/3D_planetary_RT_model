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

#look for command-line macro to 
IF RT_FLOAT:
    ctypedef float Real
    realconvert = np.float32
ELSE:
    ctypedef double Real
    realconvert = np.float64

cdef extern from "atm/chamberlain_exosphere.hpp":
    cpdef cppclass Temp_converter:
        Real lc_from_T(Real T)
        Real eff_from_T(Real T)

cdef extern from "observation_fit.hpp":
    cpdef cppclass observation_fit:
        observation_fit()
        void add_observation(vector[vector[Real]] MSO_locations, vector[vector[Real]] MSO_directions)
        void set_g_factor(Real &g)
        void generate_source_function(Real nHexo, Real Texo)
        void generate_source_function_lc(Real nHexo, Real lc)
        void generate_source_function_effv(Real nHexo, Real effv)
        vector[Real] brightness()
        Temp_converter Tconv
        
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
                loc_vec[i][j] = realconvert(loc_arr[i,j])
                dir_vec[i][j] = realconvert(dir_arr[i,j])
        self.thisptr.add_observation(loc_vec,dir_vec)
    def set_g_factor(self, Real g):
        self.thisptr.set_g_factor(g)
    def lc_from_T(self, T):
        return self.thisptr.Tconv.lc_from_T(realconvert(T))
    def eff_from_T(self, T):
        return self.thisptr.Tconv.eff_from_T(realconvert(T))
    def generate_source_function(self, Real nH, Real T):
        self.thisptr.generate_source_function(nH,T)
    def generate_source_function_lc(self, Real nH, Real lc):
        self.thisptr.generate_source_function_lc(nH,lc)
    def generate_source_function_effv(self, Real nH, Real effv):
        self.thisptr.generate_source_function_effv(nH,effv)
    def brightness(self):
        return np.asarray(self.thisptr.brightness())
