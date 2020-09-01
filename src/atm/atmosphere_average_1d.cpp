#include "atmosphere_average_1d.hpp"
using std::vector;

atmosphere_average_1d::atmosphere_average_1d() :  spherical(true), average_init(false) {}

void atmosphere_average_1d::setup() {
  //we need to integrate the relevant quantites so we can compute
  //averages

  //throughout we use n_species, n_absorber, and Temp, which are
  //defined by sibling code. This setup cannot go in the constructor
  //for this reason because the vtable will not be built and these
  //functions will not find the right implementation

  //check if the atmosphere is ready to be integrated
  if (!init)
    return;

  
  //integrate from the top of the atmosphere down to minimize floating
  //point subtraction errors
  Real log_r_int_step = (log((rmax-rMars)/r_int_scale) - log((rmin-rMars)/r_int_scale))/(n_int_steps - 1.);
  log_r_int.push_back(log((rmax-rMars)*(1-ABS)/r_int_scale));

  n_species_int.push_back(0);
  n_species_int_spherical.push_back(0);
  n_absorber_int.push_back(0);
  n_absorber_int_spherical.push_back(0);
  Tint.push_back(0);
  Tint_spherical.push_back(0);
  for (int i_int=1; i_int<n_int_steps; i_int++) {
    log_r_int.push_back(log_r_int[0]-i_int*log_r_int_step);
    if (exp(log_r_int.back()) < (rmin-rMars)/r_int_scale)
      log_r_int.back() = log((rmin-rMars)*(1+ABS)/r_int_scale);

    //scaled quantities
    Real r0s = exp(log_r_int[i_int-1])+rMars/r_int_scale;
    Real r1s = exp(log_r_int[i_int])+rMars/r_int_scale;
    Real drs = r0s-r1s;
    //unscaled quantities
    Real r0 = exp(log_r_int[i_int-1])*r_int_scale+rMars;
    Real r1 = exp(log_r_int[i_int])*r_int_scale+rMars;
    //Real dr = r0-r1;

    n_species_int.push_back( (n_species(r1) + n_species(r0))/2.0 * drs + n_species_int.back() );
    n_species_int_spherical.push_back( (n_species(r1)*r1s*r1s + n_species(r0)*r0s*r0s)/2.0 * drs + n_species_int_spherical.back() );

    n_absorber_int.push_back( (n_absorber(r1) + n_absorber(r0))/2.0 * drs + n_absorber_int.back() );
    n_absorber_int_spherical.push_back( (n_absorber(r1)*r1s*r1s + n_absorber(r0)*r0s*r0s)/2.0 * drs + n_absorber_int_spherical.back() );

    Real T0 = Temp(r0);
    Real T1 = Temp(r1);

    Tint.push_back( (T1 + T0)/2.0 * drs + Tint.back() );
    Tint_spherical.push_back( (T1*r1s*r1s + T0*r0s*r0s)/2.0 * drs + Tint_spherical.back() );
  }
  n_species_int_spline = cardinal_cubic_b_spline<Real>(n_species_int.rbegin(),
						       n_species_int.rend(),
						       log_r_int.back(),
						       log_r_int_step);
  n_species_int_spline_spherical = cardinal_cubic_b_spline<Real>(n_species_int_spherical.rbegin(),
								 n_species_int_spherical.rend(),
								 log_r_int.back(),
								 log_r_int_step);
  n_absorber_int_spline = Linear_interp(log_r_int, n_absorber_int);
  n_absorber_int_spline_spherical = Linear_interp(log_r_int, n_absorber_int_spherical);
  
  Tint_spline = Linear_interp(log_r_int, Tint);
  Tint_spline_spherical = Linear_interp(log_r_int, Tint_spherical);

  average_init=true;
}


Real atmosphere_average_1d::ravg(const Real &r0, const Real &r1,
			 const Real &q0, const Real &q1) const {
  //compute average values from integral quantities q;
  if (spherical) {
    return -3*( q1 - q0 )/( r1*r1*r1 - r0*r0*r0 );// volume is
						  // int_r0^r1(r^2)=1/3(r1^3-r0^3)
						  // (no 4pi b/c integral is r only)
    //     ^this (-) is here because integration is from the top of
    //     the atmosphere down, to minimize the chance of floating
    //     point errors
  } else {
    return -( q1 - q0 )/( r1 - r0 );
    //     ^this (-) is here because integration is from the top of
    //     the atmosphere down, to minimize the chance of floating
    //     point errors
  }
}

Real atmosphere_average_1d::n_absorber_avg(const Real &r0, const Real &r1) const {
  if (!average_init)
    assert(false && "atmosphere_average_1d not properly initialized!");
  if (spherical) {
    return ravg(r0/r_int_scale, r1/r_int_scale,
		n_absorber_int_spline_spherical(log((r0-rMars)/r_int_scale)),
		n_absorber_int_spline_spherical(log((r1-rMars)/r_int_scale)));
  } else {
    return ravg(r0/r_int_scale, r1/r_int_scale,
		n_absorber_int_spline(log((r0-rMars)/r_int_scale)),
		n_absorber_int_spline(log((r1-rMars)/r_int_scale)));
  }
}
void atmosphere_average_1d::n_absorber(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  ret_avg = n_absorber_avg(vox.rbounds[0], vox.rbounds[1]);
  assert(!isnan(ret_avg) && ret_avg >= 0 && "densities must be real and positive");
  ret_pt  = n_absorber(vox.pt.r);
  assert(!isnan(ret_pt) && ret_pt >= 0 && "densities must be real and positive");
}

Real atmosphere_average_1d::n_species_avg(const Real &r0, const Real &r1) const {
  if (!average_init)
    assert(false && "atmosphere_average_1d not properly initialized!");
  if (spherical) {
    return ravg(r0/r_int_scale, r1/r_int_scale,
		n_species_int_spline_spherical(log((r0-rMars)/r_int_scale)),
		n_species_int_spline_spherical(log((r1-rMars)/r_int_scale)));
  } else {
    return ravg(r0/r_int_scale, r1/r_int_scale,
		n_species_int_spline(log((r0-rMars)/r_int_scale)),
		n_species_int_spline(log((r1-rMars)/r_int_scale)));
  }
}
void atmosphere_average_1d::n_species(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  ret_avg = n_species_avg(vox.rbounds[0],vox.rbounds[1]);
  assert(!isnan(ret_avg) && ret_avg >= 0 && "densities must be real and positive");
  ret_pt  = n_species(vox.pt.r);
  assert(!isnan(ret_pt) && ret_pt >= 0 && "densities must be real and positive");
}


Real atmosphere_average_1d::Temp_avg(const Real &r0, const Real &r1) const {
  if (!average_init)
    assert(false && "atmosphere_average_1d not properly initialized!");
  if (spherical) {
    return ravg(r0/r_int_scale, r1/r_int_scale,
		  Tint_spline_spherical(log((r0-rMars)/r_int_scale)),
		  Tint_spline_spherical(log((r1-rMars)/r_int_scale)));
  } else {
    return ravg(r0/r_int_scale, r1/r_int_scale,
		  Tint_spline(log((r0-rMars)/r_int_scale)),
		  Tint_spline(log((r1-rMars)/r_int_scale)));
  }
}
void atmosphere_average_1d::Temp(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const {
  ret_avg = Temp_avg(vox.rbounds[0], vox.rbounds[1]);
  assert(!isnan(ret_avg) && ret_avg >= 0 && "temperatures must be real and positive");
  ret_pt = Temp(vox.pt.r);
  assert(!isnan(ret_pt) && ret_pt >= 0 && "temperatures must be real and positive");
}
