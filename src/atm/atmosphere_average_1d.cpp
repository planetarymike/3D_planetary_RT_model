#include "atmosphere_average_1d.hpp"
#include <iostream>
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


  //reset all our vectors
  log_r_int.clear();
  n_species_int.clear();
  n_species_int_spherical.clear();
  n_absorber_int.clear();
  n_absorber_int_spherical.clear();
  Tint.clear();
  Tint_spherical.clear();
  
  //integrate from the top of the atmosphere down to minimize floating
  //point subtraction errors
  double log_r_int_step = (log((rmax-rMars)/r_int_scale) - log((rmin-rMars)/r_int_scale))/(n_int_steps - 1.);
  log_r_int.push_back(log((rmax-rMars)*(1-ATMEPS)/r_int_scale));

  n_species_int.push_back(0);
  n_species_int_spherical.push_back(0);
  n_absorber_int.push_back(0);
  n_absorber_int_spherical.push_back(0);
  Tint.push_back(0);
  Tint_spherical.push_back(0);
  for (int i_int=1; i_int<n_int_steps; i_int++) {
    log_r_int.push_back(log_r_int[0]-i_int*log_r_int_step);
    if (exp(log_r_int.back()) < (rmin-rMars)/r_int_scale)
      log_r_int.back() = log((rmin-rMars)*(1+ATMEPS)/r_int_scale);

    //unscaled quantities
    double r0 = exp(log_r_int[i_int-1])*r_int_scale+rMars;
    double r1 = exp(log_r_int[i_int])*r_int_scale+rMars;
    //double dr = r0-r1;

    //scaled quantities
    double r0s = r0/r_int_scale;
    double r1s = r1/r_int_scale;
    
    add_integrated(n_species_int,
		   n_species(r0), n_species(r1),
		   r0s,           r1s,
		   /*spherical = */false);

    add_integrated(n_species_int_spherical,
		   n_species(r0), n_species(r1),
		   r0s,           r1s,
		   /*spherical = */true);

    //absorbers
    add_integrated(n_absorber_int,
		   n_absorber(r0), n_absorber(r1),
		   r0s,           r1s,
		   /*spherical = */false);

    add_integrated(n_absorber_int_spherical,
		   n_absorber(r0), n_absorber(r1),
		   r0s,           r1s,
		   /*spherical = */true);

    //temperature
    add_integrated(Tint,
		   Temp(r0), Temp(r1),
		   r0s,           r1s,
		   /*spherical = */false);

    add_integrated(Tint_spherical,
		   Temp(r0), Temp(r1),
		   r0s,           r1s,
		   /*spherical = */true);
  }
  n_species_int_spline = cardinal_cubic_b_spline<double>(n_species_int.rbegin(),
							 n_species_int.rend(),
							 log_r_int.back(),
							 log_r_int_step);
  n_species_int_spline_spherical = cardinal_cubic_b_spline<double>(n_species_int_spherical.rbegin(),
								   n_species_int_spherical.rend(),
								   log_r_int.back(),
								   log_r_int_step);
  n_absorber_int_spline = Linear_interp<double>(log_r_int, n_absorber_int);
  n_absorber_int_spline_spherical = Linear_interp<double>(log_r_int, n_absorber_int_spherical);
  
  Tint_spline = Linear_interp<double>(log_r_int, Tint);
  Tint_spline_spherical = Linear_interp<double>(log_r_int, Tint_spherical);

  average_init=true;
}

void atmosphere_average_1d::add_integrated(vector<double> &vec,
					   double q0, double q1,
					   double r0s, double r1s,
					   bool spherical) {

  double rfactor0, rfactor1;
  double diff;
  
  if (spherical) {
    rfactor0 = r0s * r0s;
    rfactor1 = r1s * r1s;
  } else {
    rfactor0 = 1.0;
    rfactor1 = 1.0;
  }

  double rfactorrel = rfactor1/rfactor0;
  double drs = r0s*(1.0-r1s/r0s);

  if (q0!=0 && q1!=0) {
    // add up relative values to avoid rounding errors
    double qrel = q1/q0;
    diff = q0*rfactor0*(1.0 + qrel*rfactorrel)/2.0 * drs;
  } else {
    // at least one of the values is zero, need to use absolute addition
    diff = rfactor0*(q0 + q1*rfactorrel)/2.0 * drs;
  }

  vec.push_back( diff + vec.back() );

  assert(!isnan(vec[vec.size()-1]) &&
	 !isnan(vec[vec.size()-2]) &&
	 !isnan(diff) &&
	 (vec[vec.size()-2]!=vec[vec.size()-1] || diff==0)
	 && "check for nans or no difference (may indicate floating point error)");
}

double atmosphere_average_1d::ravg(const double &r0, const double &r1,
				   const double &q0, const double &q1) const {
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

double atmosphere_average_1d::n_absorber_avg(const double &r0, const double &r1) const {
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

double atmosphere_average_1d::n_species_avg(const double &r0, const double &r1) const {
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


double atmosphere_average_1d::Temp_avg(const double &r0, const double &r1) const {
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
void atmosphere_average_1d::Temp(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  ret_avg = Temp_avg(vox.rbounds[0], vox.rbounds[1]);
  assert(!isnan(ret_avg) && ret_avg >= 0 && "temperatures must be real and positive");
  ret_pt = Temp(vox.pt.r);
  assert(!isnan(ret_pt) && ret_pt >= 0 && "temperatures must be real and positive");
}
