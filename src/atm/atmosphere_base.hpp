#ifndef __ATMOSPHERE_BASE_H_
#define __ATMOSPHERE_BASE_H_

#include "Real.hpp"
#include "atmo_vec.hpp"

struct atmosphere {
  doubReal rmin;// cm, minimum altitude in model atmosphere
  doubReal rexo;// cm, exobase altitude. Needs to remain unless you want
                  // to rewrite rmethod_altitude in grid_plane_parallel and
                  // grid_spherical_azimuthally_symmetric
  doubReal rmax;// cm, max altitude in model atmosphere

  bool init;//whether n_species, r_from_n_species, n_absorber, and
	    //Temp are ready to use
  
  //species density at a given radius (either subsolar or average)
  virtual doubReal n_species(const doubReal &r) const = 0;

  //returns radius (at subsolar point, or average) from species density
  virtual doubReal r_from_n_species(const doubReal &n_species) const = 0;

  //atmosphere temperature
  virtual doubReal Temp(const doubReal &r) const = 0;

  //absorber density
  virtual doubReal n_absorber(const doubReal &r) const = 0; 
  
  //function for cross sections should be defined for use with RT
  //code, but this is not required as some species (H) have multiple
  //emissions with different cross sections
  //
  //virtual doubReal species_sigma(const doubReal &T) const = 0; virtual
  //doubReal absorber_sigma(const doubReal &T) const = 0;

  atmosphere(doubReal rminn, doubReal rexoo, doubReal rmaxx);

  virtual ~atmosphere() = default;
};

#endif
