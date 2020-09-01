#ifndef __ATMOSPHERE_BASE_H_
#define __ATMOSPHERE_BASE_H_

#include "Real.hpp"
#include "atmo_vec.hpp"

struct atmosphere {
  Real rmin;// cm, minimum altitude in model atmosphere
  Real rexo;// cm, exobase altitude. Needs to remain unless you want
                  // to rewrite rmethod_altitude in grid_plane_parallel and
                  // grid_spherical_azimuthally_symmetric
  Real rmax;// cm, max altitude in model atmosphere

  bool init;//whether n_species, r_from_n_species, n_absorber, and
	    //Temp are ready to use
  
  //species density at a given radius (either subsolar or average)
  virtual Real n_species(const Real &r) const = 0;

  //returns radius (at subsolar point, or average) from species density
  virtual Real r_from_n_species(const Real &n_species) const = 0;

  //atmosphere temperature
  virtual Real Temp(const Real &n_species) const = 0;

  //absorber density
  virtual Real n_absorber(const Real &r) const = 0; 
  
  //function for cross sections should be defined for use with RT
  //code, but it not required as some species (H) have multiple
  //emissions with different cross sections
  //
  //virtual Real species_sigma(const Real &T) const = 0; virtual
  //Real absorber_sigma(const Real &T) const = 0;

  atmosphere(Real rminn, Real rexoo, Real rmaxx);
};

#endif
