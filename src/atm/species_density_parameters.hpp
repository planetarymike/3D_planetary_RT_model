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
  const double mass; // species mass
  const double alpha; // thermal diffusion coefficient
  const int n_vars; // number of variables tracked in the thermosphere
  diffusion_coefs diff; // object to calculate diffusion coefficients

  double escape_flux; // boundary condition for differential equation
  double rexo;
  temperature *temp; // pointer to temperature, needed by operator () to get derivatives

  species_density_parameters(const double masss,
			     const double alphaa,
			     const int n_varss,
			     diffusion_coefs difff);


  // used by get_thermosphere_densities()
  virtual void operator()( const vector<double> &x , vector<double> &dxdr , const double &r ) = 0;
  // returns the diffeq derivatives of the thermosphere diffusion
  // equation for log(nCO2) and the species of interest
  // note: log(nCO2) must always be the zeroth element of x
  //       log(n_species) must always be the last element of x

  // called by thermosphere_exosphere to populate thermosphere interpolation arrays
  virtual void get_thermosphere_density_arrays(const double &n_species_exo,
					       const double &n_CO2_exo,
					       const double &rexo,
					       const double &rmin,

					       const double &escape_fluxx,
					       temperature *tempp,

					       vector<double> &log_nCO2_thermosphere,
					       vector<double> &log_n_species_thermosphere,
					       vector<double> &r_thermosphere,
					       const int n_thermosphere_steps,
					       bool get_interpolation_points) = 0;

  //routines to get exact densities
  double get_CO2_exobase_density(const double &nCO2_lower,
				 const double &r_lower,
				 const double &rexoo,
				 temperature *tempp);

  void get_thermosphere_density_exact(const double &n_species_exo,
				      const double &n_CO2_exo,
				      const double &rexoo,
				      
				      const double &escape_fluxx,
				      temperature *tempp,
				      
				      const double &r, 
				      double &nCO2,
				      double &n_species);
};


struct hydrogen_density_parameters : species_density_parameters {
  static constexpr double DH0_hydrogen = 8.4e17;// cm^2 s^-1
  static constexpr double s_hydrogen = 0.6;
  static constexpr double alpha_hydrogen = -0.25; // thermal diffusion coefficient

  hydrogen_density_parameters();
  
  // returns the diffusion equation derivatives
  void operator()( const vector<double> &x , vector<double> &dxdr , const double &r ) override;

  // called by thermosphere_exosphere to populate thermosphere interpolation arrays
  void get_thermosphere_density_arrays(const double &n_species_exo,
				       const double &n_CO2_exo,
				       const double &rexoo,
				       const double &rmin,
				       
				       const double &escape_fluxx,
				       temperature *tempp,
				       
				       vector<double> &log_nCO2_thermosphere,
				       vector<double> &log_n_species_thermosphere,
				       vector<double> &r_thermosphere,
				       const int n_thermosphere_steps,
				       bool get_interpolation_points) override;
};

#endif
