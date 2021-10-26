#ifndef __SPECIES_DENSITY_PARAMETERS_H
#define __SPECIES_DENSITY_PARAMETERS_H

#include "temperature.hpp"
#include "diffusion.hpp"

#include <vector>
using std::vector;

// define species properties, including mass, diffusion coefficients,
// and the differential equation used to solve for the CO2 and species
// number densities in the themosphere

struct species_density_parameters {
  const doubReal mass; // species mass
  const doubReal alpha; // thermal diffusion coefficient
  diffusion_coefs diff; // object to calculate diffusion coefficients

  doubReal escape_flux; // boundary condition for differential equation
  doubReal rexo;
  temperature *temp; // pointer to temperature, needed by operator () to get derivatives

  species_density_parameters(const doubReal masss,
			     const doubReal alphaa,
			     diffusion_coefs difff);


  // used by get_thermosphere_densities()
  virtual void operator()( const vector<doubReal> &x , vector<doubReal> &dxdr , const doubReal &r ) = 0;
  // returns the diffeq derivatives of the thermosphere diffusion
  // equation for log(nCO2) and the species of interest
  // note: log(nCO2) must always be the zeroth element of x
  //       log(n_species) must always be the last element of x

  // called by thermosphere_exosphere to populate thermosphere interpolation arrays
  virtual void get_thermosphere_density_arrays(const doubReal &n_species_exo,
					       const doubReal &n_CO2_exo,
					       const doubReal &rexo,
					       const doubReal &rmin,

					       const doubReal &escape_fluxx,
					       temperature *tempp,

					       vector<doubReal> &log_nCO2_thermosphere,
					       vector<doubReal> &log_n_species_thermosphere,
					       vector<doubReal> &r_thermosphere,
					       const int n_thermosphere_steps,
					       bool get_interpolation_points) = 0;

  //routines to get exact densities
  doubReal get_CO2_exobase_density(const doubReal &nCO2_lower,
				 const doubReal &r_lower,
				 const doubReal &rexoo,
				 temperature *tempp);

  void get_thermosphere_density_exact(const doubReal &n_species_exo,
				      const doubReal &n_CO2_exo,
				      const doubReal &rexoo,
				      
				      const doubReal &escape_fluxx,
				      temperature *tempp,
				      
				      const doubReal &r, 
				      doubReal &nCO2,
				      doubReal &n_species);
};


struct hydrogen_density_parameters : species_density_parameters {
  static constexpr doubReal DH0_hydrogen = 8.4e17;// cm^2 s^-1
  static constexpr doubReal s_hydrogen = 0.6;
  static constexpr doubReal alpha_hydrogen = -0.25; // thermal diffusion coefficient

  hydrogen_density_parameters();
  
  // returns the diffusion equation derivatives
  void operator()( const vector<doubReal> &x , vector<doubReal> &dxdr , const doubReal &r ) override;

  // called by thermosphere_exosphere to populate thermosphere interpolation arrays
  void get_thermosphere_density_arrays(const doubReal &n_species_exo,
				       const doubReal &n_CO2_exo,
				       const doubReal &rexoo,
				       const doubReal &rmin,
				       
				       const doubReal &escape_fluxx,
				       temperature *tempp,
				       
				       vector<doubReal> &log_nCO2_thermosphere,
				       vector<doubReal> &log_n_species_thermosphere,
				       vector<doubReal> &r_thermosphere,
				       const int n_thermosphere_steps,
				       bool get_interpolation_points) override;
};

struct oxygen_density_parameters : species_density_parameters {
  static constexpr doubReal DH0_oxygen = 4.4e17;// cm^2 s^-1
  static constexpr doubReal s_oxygen = 0.5;
  static constexpr doubReal alpha_oxygen = 0.0; // thermal diffusion coefficient

  oxygen_density_parameters();
  
  // returns the diffusion equation derivatives
  void operator()( const vector<doubReal> &x , vector<doubReal> &dxdr , const doubReal &r ) override;

  // called by thermosphere_exosphere to populate thermosphere interpolation arrays
  void get_thermosphere_density_arrays(const doubReal &n_species_exo,
				       const doubReal &n_CO2_exo,
				       const doubReal &rexoo,
				       const doubReal &rmin,
				       
				       const doubReal &escape_fluxx,
				       temperature *tempp,
				       
				       vector<doubReal> &log_nCO2_thermosphere,
				       vector<doubReal> &log_n_species_thermosphere,
				       vector<doubReal> &r_thermosphere,
				       const int n_thermosphere_steps,
				       bool get_interpolation_points) override;
};


#endif
