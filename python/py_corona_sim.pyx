# py_corona_sim.pyx -- Cython wrapper for H corona simulation object
# distutils: language = c++
# cython: language_level=3
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
import importlib
import os

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

cdef extern from "chamberlain_exosphere.hpp":
    cdef cppclass Temp_converter:
        Real lc_from_T(Real T)
        Real eff_from_T(Real T)

# name of Quemerais IPH source function, packaged along with *.so file
iph_sfn_basename = 'quemerais_IPH_sourcefn_fsm99td12v20t80.dat' # basename of source file
# fully qualified name determined when class Pyobservation_fit is created
        
cdef extern from "observation_fit.hpp":
    cdef cppclass observation_fit:
        observation_fit(string sfn_fname)

        void add_observation(vector[vector[Real]] MSO_locations, vector[vector[Real]] MSO_directions)

        void set_g_factor(vector[Real] &g)

        void simulate_iph(bool sim_iphh);
        void add_observation_ra_dec(vector[Real] mars_ecliptic_coords,
			            vector[Real] &RAA,
	                            vector[Real] &Decc);

        void generate_source_function(Real nHexo, Real Texo,
                                      string atmosphere_fname,
                                      string sourcefn_fname,
                                      bool plane_parallel,
                                      bool deuterium)
        void generate_source_function_lc(Real nHexo, Real lc, 
                                      string atmosphere_fname,
                                      string sourcefn_fname,
                                      bool plane_parallel,
                                      bool deuterium)
        void generate_source_function_effv(Real nHexo, Real effv,
                                           string atmosphere_fname,
                                           string sourcefn_fname,
                                           bool plane_parallel,
                                           bool deuterium)

        void generate_source_function_variable_thermosphere(Real nHexo, Real Texo,
                                                            Real nCO2rminn,
						            Real rexoo,
						            Real rminn,
						            Real rmaxx,
						            Real rmindiffusionn,
						            #extra args for krasnopolsky_temp					  
						            Real T_tropo,
						            Real r_tropo,
						            Real shape_parameter,				  
                                                            string atmosphere_fname,
                                                            string sourcefn_fname,
                                                            bool plane_parallel,
                                                            bool deuterium)
        
        void generate_source_function_nH_asym(Real nHexo, Real Texo,
                                              Real asym,
                                              string sourcefn_fname,
                                              bool deuterium)

        void generate_source_function_temp_asym(Real nHavg,
                                                Real Tnoon, Real Tmidnight,
                                                string sourcefn_fname,
                                                bool deuterium)
        
        void generate_source_function_temp_asym(Real nHavg,
                                                Real Tnoon, Real Tmidnight,
						Real nCO2rminn,
						Real rexoo,
						Real rminn,
						Real rmaxx,
						Real rmindiffusionn,
						#extra args for krasnopolsky_temp					  
						Real T_tropo,
						Real r_tropo,
						Real shape_parameter,
                                                #power for temperature in the expression n*T^p = const.
                                                Real Tpower,
                                                string sourcefn_fname,
                                                bool deuterium)

        void generate_source_function_tabular_atmosphere(Real rmin, Real rexo, Real rmax,
							 vector[double] &alt_nH, vector[double] &log_nH,
							 vector[double] &alt_nCO2, vector[double] &log_nCO2,
							 vector[double] &alt_temp, vector[double] &temp,
							 bool compute_exosphere,
                                                         bool plane_parallel,
                                                         bool deuterium,
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

        Real get_CO2_exobase_density()
        void reset_CO2_exobase_density()
        void set_CO2_exobase_density(Real nCO2)

        void save_influence_matrix(string fname)
        void save_influence_matrix_O_1026(string fname)
        
        vector[vector[Real]] brightness()
        vector[vector[Real]] species_col_dens()
        vector[vector[Real]] tau_species_final()
        vector[vector[Real]] tau_absorber_final()
        vector[vector[Real]] iph_brightness_observed()
        vector[vector[Real]] iph_brightness_unextincted()
        Temp_converter Tconv

        vector[vector[Real]] D_brightness()
        vector[vector[Real]] D_col_dens()
        vector[vector[Real]] tau_D_final()

        void O_1026_generate_source_function(Real nOexo,
                                             Real Texo,
                                             Real solar_brightness_lyman_beta, # 1.69e-3 is a good number for solar minimum
                                             string atmosphere_fname,
                                             string sourcefn_fname)
        vector[vector[Real]] O_1026_brightness()

        void lyman_multiplet_generate_source_function(Real nHexo,
                                                      Real Texo,
                                                      # Real solar_brightness_lyman_alpha,
                                                      # Real solar_brightness_lyman_beta,
                                                      string atmosphere_fname,
                                                      string sourcefn_fname)
        vector[vector[Real]] lyman_multiplet_brightness()


        
cdef class Pyobservation_fit:
    cdef observation_fit *thisptr #holds the reference to the cpp class
    def __cinit__(self):
        # if we've gotten here, the module must be known to python,
        # find its origin file
        module_spec_gpu = importlib.util.find_spec('py_corona_sim_gpu')
        module_spec_cpu = importlib.util.find_spec('py_corona_sim_cpu')
        if module_spec_gpu is not None:
            module_fname = module_spec_gpu.origin
        elif module_spec_cpu is not None:
            module_fname = module_spec_cpu.origin
        else:
            raise RuntimeError("Cannot locate py_corona_sim package source file")

        module_dirname = os.path.dirname(module_fname)
        iph_sfn_fname = os.path.join(module_dirname, iph_sfn_basename)

        # check if the file exists
        if not os.path.exists(iph_sfn_fname):
            raise FileNotFoundError("Cannot find necessary "
                                    +"IPH source function file "
                                    +f"{iph_sfn_fname}, check "
                                    +"that it exists.")

        print("loading IPH source function from file " + iph_sfn_fname)
        self.thisptr = new observation_fit(iph_sfn_fname.encode('utf-8'))
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

    def simulate_iph(self,
                     bool sim_iphh):
        self.thisptr.simulate_iph(sim_iphh)

        
    def add_observation_ra_dec(self,
                               mars_ecliptic_coords_arr,
			       RA_arr,
	                       Dec_arr):

        cdef vector[Real] marspos, RA, Dec

        marspos.resize(3)
        for i in range(3):
            marspos[i] = realconvert(mars_ecliptic_coords_arr[i])
        
        RA.resize(RA_arr.shape[0])
        Dec.resize(Dec_arr.shape[0])
        for i in range(RA_arr.shape[0]):
            RA[i]  = realconvert(RA_arr[i])
            Dec[i] = realconvert(Dec_arr[i])
        
        self.thisptr.add_observation_ra_dec(marspos,
                                            RA,
                                            Dec)
        
    def lc_from_T(self, T):
        return self.thisptr.Tconv.lc_from_T(realconvert(T))
    def eff_from_T(self, T):
        return self.thisptr.Tconv.eff_from_T(realconvert(T))

    def generate_source_function(self, Real nH, Real T,
                                 atmosphere_fname = "",
                                 sourcefn_fname = "",
                                 plane_parallel = False,
                                 deuterium = False):
        self.thisptr.generate_source_function(nH,T,
                                              atmosphere_fname.encode('utf-8'),
                                              sourcefn_fname.encode('utf-8'),
                                              plane_parallel,
                                              deuterium)
    def generate_source_function_lc(self, Real nH, Real lc,
                                    atmosphere_fname = "",
                                    sourcefn_fname = "",
                                    plane_parallel = False,
                                    deuterium = False):
        self.thisptr.generate_source_function_lc(nH,lc,
                                                 atmosphere_fname.encode('utf-8'),
                                                 sourcefn_fname.encode('utf-8'),
                                                 plane_parallel,
                                                 deuterium)
    def generate_source_function_effv(self, Real nH, Real effv,
                                      atmosphere_fname = "",
                                      sourcefn_fname = "",
                                      plane_parallel = False,
                                      deuterium = False):
        self.thisptr.generate_source_function_effv(nH,effv,
                                                   atmosphere_fname.encode('utf-8'),
                                                   sourcefn_fname.encode('utf-8'),
                                                   plane_parallel,
                                                   deuterium)
    def generate_source_function_variable_thermosphere(self,
                                                       Real nHexo, Real Texo,
                                                       Real nCO2rminn,
						       Real rexoo,
						       Real rminn,
						       Real rmaxx,
						       Real rmindiffusionn,
						       #extra args for krasnopolsky_temp					  
						       Real T_tropo,
						       Real r_tropo,
						       Real shape_parameter,				  
                                                       atmosphere_fname = "",
                                                       sourcefn_fname = "",
                                                       plane_parallel = False,
                                                       deuterium = False):
        self.thisptr.generate_source_function_variable_thermosphere(nHexo, Texo,
                                                                    nCO2rminn,
						                    rexoo,
						                    rminn,
						                    rmaxx,
						                    rmindiffusionn,
						                    #extra args for krasnopolsky_temp					  
						                    T_tropo,
						                    r_tropo,
						                    shape_parameter,
				                                    atmosphere_fname.encode('utf-8'),
                                                                    sourcefn_fname.encode('utf-8'),
                                                                    plane_parallel,
                                                                    deuterium)
        
    def generate_source_function_nH_asym(self, Real nH, Real Texo,
                                         Real asym,
                                         sourcefn_fname = "",
                                         deuterium = False):
        self.thisptr.generate_source_function_nH_asym(nH,Texo,
                                                      asym,
                                                      sourcefn_fname.encode('utf-8'),
                                                      deuterium)

    def generate_source_function_temp_asym(self, Real nHavg,
                                           Real Tnoon, Real Tmidnight,
                                           sourcefn_fname = "",
                                           deuterium = False):
        self.thisptr.generate_source_function_temp_asym(nHavg,
                                                        Tnoon,Tmidnight,
                                                        sourcefn_fname.encode('utf-8'),
                                                        deuterium)

    def generate_source_function_temp_asym_full(self, Real nHavg,
                                                Real Tnoon, Real Tmidnight,
                                                Real nCO2rminn,
					        Real rexoo,
					        Real rminn,
					        Real rmaxx,
					        Real rmindiffusionn,
					        #extra args for krasnopolsky_temp					  
					        Real T_tropo,
					        Real r_tropo,
					        Real shape_parameter,
			                        #power for temperature in the expression n*T^p = const.
                                                Real Tpower,
                                                sourcefn_fname = "",
                                                deuterium = False):
        self.thisptr.generate_source_function_temp_asym(nHavg,
                                                        Tnoon,Tmidnight,
                                                        nCO2rminn,
					                rexoo,
					                rminn,
					                rmaxx,
					                rmindiffusionn,
					                #extra args for krasnopolsky_temp					  
					                T_tropo,
					                r_tropo,
                                                        shape_parameter,
						        #power for temperature in the expression n*T^p = const.
                                                        Tpower,
                                                        sourcefn_fname.encode('utf-8'),
                                                        deuterium)

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
                                                    deuterium = False,
                                                    sourcefn_fname = ""):
        cdef vector[double] alt_nH, log_nH, alt_nCO2, log_nCO2, alt_Temp, Temp
        alt_nH.resize(atm_dict['alt_nH'].shape[0])
        log_nH.resize(atm_dict['log_nH'].shape[0])
        for i in range(atm_dict['alt_nH'].shape[0]):
            alt_nH[i] =  np.float64(atm_dict['alt_nH'][i])
            log_nH[i] =  np.float64(atm_dict['log_nH'][i])
        alt_nCO2.resize(atm_dict['alt_nCO2'].shape[0])
        log_nCO2.resize(atm_dict['log_nCO2'].shape[0])
        for i in range(atm_dict['alt_nCO2'].shape[0]):
            alt_nCO2[i] =  np.float64(atm_dict['alt_nCO2'][i])
            log_nCO2[i] =  np.float64(atm_dict['log_nCO2'][i])
        alt_Temp.resize(atm_dict['alt_Temp'].shape[0])
        Temp.resize(atm_dict['Temp'].shape[0])
        for i in range(atm_dict['alt_Temp'].shape[0]):
            alt_Temp[i] =  np.float64(atm_dict['alt_Temp'][i])
            Temp[i] =  np.float64(atm_dict['Temp'][i])
        
        self.thisptr.generate_source_function_tabular_atmosphere(atm_dict['rmin'],
                                                                 atm_dict['rexo'],
                                                                 atm_dict['rmax'],
                                                                 alt_nH,   log_nH,
                                                                 alt_nCO2, log_nCO2,
                                                                 alt_Temp, Temp,
                                                                 compute_exosphere,
                                                                 plane_parallel,
                                                                 deuterium,
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

    def get_CO2_exobase_density(self):
        return self.thisptr.get_CO2_exobase_density();

    def reset_CO2_exobase_density(self):
        self.thisptr.reset_CO2_exobase_density();

    def set_CO2_exobase_density(self, nCO2):
        self.thisptr.set_CO2_exobase_density(nCO2);
            
    def brightness(self):
        return np.asarray(self.thisptr.brightness())

    def species_col_dens(self):
        return np.asarray(self.thisptr.species_col_dens())

    def tau_species_final(self):
        return np.asarray(self.thisptr.tau_species_final())

    def tau_absorber_final(self):
        return np.asarray(self.thisptr.tau_absorber_final())

    def iph_brightness_observed(self):
        return np.asarray(self.thisptr.iph_brightness_observed())

    def iph_brightness_unextincted(self):
        return np.asarray(self.thisptr.iph_brightness_unextincted())

    def D_brightness(self):
        return np.asarray(self.thisptr.D_brightness())

    def D_col_dens(self):
        return np.asarray(self.thisptr.D_col_dens())

    def tau_D_final(self):
        return np.asarray(self.thisptr.tau_D_final())

    def save_influence_matrix(self, fname):
        self.thisptr.save_influence_matrix(fname.encode('utf-8'))

    def save_influence_matrix_O_1026(self, fname):
        self.thisptr.save_influence_matrix_O_1026(fname.encode('utf-8'))

    
    def O_1026_generate_source_function(self,
                                        Real nO,
                                        Real T,
                                        Real solar_brightness_lyman_beta,
                                        atmosphere_fname = "",
                                        sourcefn_fname = ""):
        self.thisptr.O_1026_generate_source_function(nO,
                                                     T,
                                                     solar_brightness_lyman_beta,
                                                     atmosphere_fname.encode('utf-8'),
                                                     sourcefn_fname.encode('utf-8'))
    def O_1026_brightness(self):
        return np.asarray(self.thisptr.O_1026_brightness())

    def lyman_multiplet_generate_source_function(self,
                                                 Real nH,
                                                 Real T,
                                                 # Real solar_brightness_lyman_alpha,
                                                 # Real solar_brightness_lyman_beta,
                                                 atmosphere_fname = "",
                                                 sourcefn_fname = ""):
        self.thisptr.lyman_multiplet_generate_source_function(nH,
                                                              T,
                                                              # solar_brightness_lyman_alpha,
                                                              # solar_brightness_lyman_beta,
                                                              atmosphere_fname.encode('utf-8'),
                                                              sourcefn_fname.encode('utf-8'))

    def lyman_multiplet_brightness(self):
        return np.asarray(self.thisptr.lyman_multiplet_brightness())

