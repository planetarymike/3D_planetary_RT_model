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

#look for command-line macro to determine whether to use 32- or 64 bit floats
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

        void set_g_factor(vector[Real] &g)

        void generate_source_function(Real nHexo, Real Texo,
                                      string atmosphere_fname,
                                      string sourcefn_fname,
                                      bool plane_parallel)
        void generate_source_function_lc(Real nHexo, Real lc, 
                                      string atmosphere_fname,
                                      string sourcefn_fname,
                                      bool plane_parallel)
        void generate_source_function_effv(Real nHexo, Real effv,
                                           string atmosphere_fname,
                                           string sourcefn_fname,
                                           bool plane_parallel)
    
        void generate_source_function_nH_asym(Real nHexo, Real Texo,
                                              Real asym,
                                              string sourcefn_fname)

        void generate_source_function_temp_asym(Real nHavg,
                                                Real Tnoon, Real Tmidnight,
                                                string sourcefn_fname)

        void generate_source_function_tabular_atmosphere(Real rmin, Real rexo, Real rmax,
							 vector[Real] &alt_nH, vector[Real] &log_nH,
							 vector[Real] &alt_nCO2, vector[Real] &log_nCO2,
							 vector[Real] &alt_temp, vector[Real] &temp,
							 bool compute_exosphere,
                                                         bool plane_parallel,
                                                         string sourcefn_fname)

        void set_use_CO2_absorption(bool use_CO2_absorption)
        void set_use_temp_dependent_sH(bool use_temp_dependent_sH, Real constant_temp_sH)

        void set_sza_method_uniform()
        void set_sza_method_uniform_cos()
                                           
        void reset_H_lya_xsec_coef(Real xsec_coef)
        void reset_H_lya_xsec_coef() # uses C++ default
        void reset_H_lyb_xsec_coef(Real xsec_coef)
        void reset_H_lyb_xsec_coef() # uses C++ default
        void reset_CO2_lya_xsec(Real xsec)
        void reset_CO2_lya_xsec() # uses C++ default
        void reset_CO2_lyb_xsec(Real xsec)
        void reset_CO2_lyb_xsec() # uses C++ default
        
        vector[vector[Real]] brightness()
        vector[vector[Real]] tau_species_final()
        vector[vector[Real]] tau_absorber_final()
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

    def set_g_factor(self, vector[Real] g):
        self.thisptr.set_g_factor(g)

    def lc_from_T(self, T):
        return self.thisptr.Tconv.lc_from_T(realconvert(T))
    def eff_from_T(self, T):
        return self.thisptr.Tconv.eff_from_T(realconvert(T))

    def generate_source_function(self, Real nH, Real T,
                                 atmosphere_fname = "",
                                 sourcefn_fname = "",
                                 plane_parallel = False):
        self.thisptr.generate_source_function(nH,T,
                                              atmosphere_fname.encode('utf-8'),
                                              sourcefn_fname.encode('utf-8'),
                                              plane_parallel)
    def generate_source_function_lc(self, Real nH, Real lc,
                                    atmosphere_fname = "",
                                    sourcefn_fname = "",
                                    plane_parallel = False):
        self.thisptr.generate_source_function_lc(nH,lc,
                                                 atmosphere_fname.encode('utf-8'),
                                                 sourcefn_fname.encode('utf-8'),
                                                 plane_parallel)
    def generate_source_function_effv(self, Real nH, Real effv,
                                      atmosphere_fname = "",
                                      sourcefn_fname = "",
                                      plane_parallel = False):
        self.thisptr.generate_source_function_effv(nH,effv,
                                                   atmosphere_fname.encode('utf-8'),
                                                   sourcefn_fname.encode('utf-8'),
                                                   plane_parallel)

    def generate_source_function_nH_asym(self, Real nH, Real Texo,
                                      Real asym,
                                      sourcefn_fname = ""):
        self.thisptr.generate_source_function_nH_asym(nH,Texo,
                                                      asym,
                                                      sourcefn_fname.encode('utf-8'))

    def generate_source_function_temp_asym(self, Real nHavg,
                                           Real Tnoon, Real Tmidnight,
                                           sourcefn_fname = ""):
        self.thisptr.generate_source_function_temp_asym(nHavg,
                                                        Tnoon,Tmidnight,
                                                        sourcefn_fname.encode('utf-8'))

    def get_example_tabular_atmosphere(self):
        rmin = 3395e5 +    80e5
        rexo = 3395e5 +   200e5
        rmax = 3395e5 + 50000e5
        alt_example = np.linspace(rmin,rmax,10)
        
        return {'rmin':rmin,
                'rexo':rexo,
                'rmax':rmax,
                'alt_nH':alt_example,
                'log_nH':np.ones_like(alt_example),
                'alt_nCO2':alt_example,
                'log_nCO2':np.ones_like(alt_example),
                'alt_Temp':alt_example,
                'Temp':np.ones_like(alt_example)}
    
    def generate_source_function_tabular_atmosphere(self, atm_dict,
                                                    compute_exosphere = False,
                                                    plane_parallel = False,
                                                    sourcefn_fname = ""):
        cdef vector[Real] alt_nH, log_nH, alt_nCO2, log_nCO2, alt_Temp, Temp
        alt_nH.resize(atm_dict['alt_nH'].shape[0])
        log_nH.resize(atm_dict['log_nH'].shape[0])
        for i in range(atm_dict['alt_nH'].shape[0]):
            alt_nH[i] = realconvert(atm_dict['alt_nH'][i])
            log_nH[i] = realconvert(atm_dict['log_nH'][i])
        alt_nCO2.resize(atm_dict['alt_nCO2'].shape[0])
        log_nCO2.resize(atm_dict['log_nCO2'].shape[0])
        for i in range(atm_dict['alt_nCO2'].shape[0]):
            alt_nCO2[i] = realconvert(atm_dict['alt_nCO2'][i])
            log_nCO2[i] = realconvert(atm_dict['log_nCO2'][i])
        alt_Temp.resize(atm_dict['alt_Temp'].shape[0])
        Temp.resize(atm_dict['Temp'].shape[0])
        for i in range(atm_dict['alt_Temp'].shape[0]):
            alt_Temp[i] = realconvert(atm_dict['alt_Temp'][i])
            Temp[i] = realconvert(atm_dict['Temp'][i])
        
        self.thisptr.generate_source_function_tabular_atmosphere(atm_dict['rmin'],
                                                                 atm_dict['rexo'],
                                                                 atm_dict['rmax'],
                                                                 alt_nH,   log_nH,
                                                                 alt_nCO2, log_nCO2,
                                                                 alt_Temp, Temp,
                                                                 compute_exosphere,
                                                                 plane_parallel,
                                                                 sourcefn_fname.encode('utf-8'))

    def set_use_CO2_absorption(self, use_CO2_absorption = True):
        self.thisptr.set_use_CO2_absorption(use_CO2_absorption)
        
    def set_use_temp_dependent_sH(self, use_temp_dependent_sH = True, constant_temp_sH = -1):
        self.thisptr.set_use_temp_dependent_sH(use_temp_dependent_sH,constant_temp_sH)

    def set_sza_method_uniform(self):
        self.thisptr.set_sza_method_uniform()

    def set_sza_method_uniform_cos(self):
        self.thisptr.set_sza_method_uniform_cos()
        
    def reset_H_lya_xsec_coef(self, xsec_coef = None):
        if xsec_coef==None:
            #use C++ defaults
            self.thisptr.reset_H_lya_xsec_coef()
        else:
            self.thisptr.reset_H_lya_xsec_coef(xsec_coef)
    def reset_H_lyb_xsec_coef(self, xsec_coef = None):
        if xsec_coef==None:
            #use C++ defaults
            self.thisptr.reset_H_lyb_xsec_coef()
        else:
            self.thisptr.reset_H_lyb_xsec_coef(xsec_coef)
    def reset_CO2_lya_xsec(self, xsec = None):
        if xsec==None:
            #use C++ defaults
            self.thisptr.reset_CO2_lya_xsec()
        else:
            self.thisptr.reset_CO2_lya_xsec(xsec)
    def reset_CO2_lyb_xsec(self, xsec = None):
        if xsec==None:
            #use C++ defaults
            self.thisptr.reset_CO2_lyb_xsec() 
        else:
            self.thisptr.reset_CO2_lyb_xsec(xsec)

    def brightness(self):
        return np.asarray(self.thisptr.brightness())

    def tau_species_final(self):
        return np.asarray(self.thisptr.tau_species_final())

    def tau_absorber_final(self):
        return np.asarray(self.thisptr.tau_absorber_final())
