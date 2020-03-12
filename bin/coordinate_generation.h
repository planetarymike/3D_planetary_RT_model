//coordinate_generation.h -- routines to get point and ray coordinates

#ifndef __coordinate_generation_h
#define __coordinate_generation_h

#include <vector>
#include <cmath>
#include "constants.h"
#include "atmosphere.h"
#include "interp.h"
#include "gauss_legendre_quadrature.h" // numerical recipies gauss quadrature points.
// would switch to boost, but gauss weights are templated, not dynamic



using std::vector;
using std::exp;
using std::log;

void gauss_quadrature_points(vector<double> &pts,
			     vector<double> &wts,
			     double start,
			     double end,
			     int npts) {

  pts.resize(npts);
  wts.resize(npts);
  
  gauleg(start,end,pts,wts);
  
}

void uniform_quadrature_points(vector<double> &pts,
			       vector<double> &wts,
			       double start,
			       double end,
			       int npts,
			       bool cyclic = false,
			       double offset = 0.0) {

  pts.resize(npts);
  wts.resize(npts);
  
  int ndivisions;
  if (cyclic) 
    ndivisions = npts+1;
  else 
    ndivisions = npts;

  double step = ( end - start )/double(ndivisions-1);

  for (int i = 0; i<npts; i++) {
    pts[i] = start + i * step;
    if (cyclic)
      pts[i]+=step*offset;
    wts[i] = step;
  }
  
}

void get_radial_log_linear_points(vector<double> &rpts,
				  int nrpts,
				  double rminatm,
				  double rexo,
				  double rmaxatm) {
  // gets radial points that are split, with half linearly spaced below
  //   the exobase and half logarithmically spaced above. 
  
  int nbelowrexo=nrpts/2;
  
  double logmax = log(rmaxatm-rMars);
  double logmin = log(rexo-rMars);
  double logspace = (logmax-logmin)/((double) nrpts-nbelowrexo);

  double linspace = (rexo-rminatm)/((double) nbelowrexo-1);

  rpts.clear();
  for (int i = 0; i<nrpts; i++) {
    if (i<nbelowrexo) 
      rpts.push_back(rminatm+i*linspace);
    else
      rpts.push_back(
		     exp(
			 logmin + (i-nbelowrexo+1)*logspace )
		     + rMars);
  }

}




// struct atmosphere_radial_boundaries {
//   atmosphere &atm;

//   atmosphere_radial_boundaries(atmosphere & atmm) : atm(atmm) { }
  
//   vector<double> lognH(int n_bounds) {
//     using std::log;
//     using std::exp;
//     int nrsteps = 1e3;
//     double logrmin = log(atm.rmin);
//     double logrmax = log(atm.rmax);
//     double logr_step = (logrmax-logrmin)/(n_bounds-1.);

//     vector<double> logr;
//     vector<double> lognH;
//     for (int i=0;i<nrsteps;i++) {
//       logr.push_back(logrmin+i*logr_step);
//       lognH.push_back(log(atm.nH(exp(logr[i]))));
//     }

//     Linear_interp 

      
//       vlognH.push_back(log(atm.nH(r = r + km_step)));

//     cardinal_cubic_b_spline<double> lognH_spline(vlognH.rbegin(),
// 						 vlognH.rend(),
// 						 r-km_step,
// 						 -km_step);
    
//     double lognH_max = log(atm.nH(atm.rmin));
//     double lognH_min = log(atm.nH(atm.rmax));
//     double lognH_step = (lognH_max-lognH_min)/(n_bounds-1.);

//     vector<double> bounds;
//     bounds.push_back(atm.rmin);
//     for(int i=1;i<n_bounds-1;i++) {
//       double lognH_target=lognH_max-i*lognH_step;
//       auto nH_n = [=](double r) { return std::make_pair(lognH_spline(r) - lognH_target, lognH_spline.prime(r)); };
//       std::uintmax_t iterations = 1000;
//       double rtarget = boost::math::tools::newton_raphson_iterate(nH_n,
// 								  1.001*bounds[i-1], atm.rmin, atm.rmax,
// 								  20, iterations);
//       bounds.push_back(rtarget);
//     }
//     bounds.push_back(atm.rmax);

//     return bounds;

//   }

  // vector<double> colH;
  // vector<double> coltauH;
  // vector<double> colCO2;
  // vector<double> coltauCO2;
  // vector<double> coltaueff;
  // vector<double> col_radii;

  // struct column_integrator {
  //   atmosphere &atm;
    
  //   column_integrator(atmosphere &atmm) : atm(atmm) { }
    
  //   void operator()( const vector<double> &x , vector<double> &dxdr , const double &r ) {
  //     dxdr[0] = -1*atm.nH(r);
  //     dxdr[1] = dxdr[0]*atm.sH_lya(r);
  //     dxdr[2] = -1*atm.nCO2(r);
  //     dxdr[3] = dxdr[2]*atm.sCO2_lya(r);
  //     dxdr[4] = dxdr[1]*exp(-x[3]);
  //   }
  // };
  

  // void compute_column_density() {
  //   column_integrator n_int(atm);
    
  //   vector<double> col(5,0);
  //   runge_kutta4< vector<double> > stepper;

  //   colH.clear();
  //   coltauH.clear();
  //   colCO2.clear();
  //   coltauCO2.clear();
  //   coltaueff.clear();
  //   col_radii.clear();
    
  //   integrate_const( stepper , n_int,
  // 		     col , atm.rmax , atm.rmin , km_step,
  // 		     push_back_quantities( &colH,
  // 					   &coltauH,
  // 					   &colCO2,
  // 					   &coltauCO2,
  // 					   &coltaueff,
  // 					   &col_radii)
  // 		     );
  // }


  // vector<double> tauH(int n_bounds) {
  //   compute_column_density();

  //   cardinal_cubic_b_spline<double> tauH_spline(coltauH.rbegin(),
  // 						coltauH.rend(),
  // 						col_radii.back(),
  // 						-km_step);
    
  //   double tauH_min = 0;
  //   double tauH_max = coltauH.back();
  //   double tauH_step = (tauH_max-tauH_min)/(n_bounds-1.);

  //   vector<double> bounds;
  //   bounds.push_back(atm.rmin);
  //   for (int i=1;i<n_bounds-1;i++) {
  //     double tau_target = tauH_max - i*tauH_step;
  //     auto tau_n = [=](double r) { return std::make_pair(tauH_spline(r) - tau_target, tauH_spline.prime(r)); };
  //     std::uintmax_t iterations = 1000;
  //     double rtarget = boost::math::tools::newton_raphson_iterate(tau_n,
  // 								  1.01*bounds[i-1], atm.rmin, atm.rmax,
  // 								  20, iterations);
  //     bounds.push_back(rtarget);
  //   }
  //   bounds.push_back(atm.rmax);

  //   return bounds;
  // }

  
  // vector<double> taueff(int n_bounds) {
  //   compute_column_density();

  //   cardinal_cubic_b_spline<double> taueff_spline(coltaueff.rbegin(),
  // 						  coltaueff.rend(),
  // 						  col_radii.back(),
  // 						  -km_step);
    
  //   double taueff_min = 0;
  //   double taueff_max = coltaueff.back();
  //   double taueff_step = (taueff_max-taueff_min)/(n_bounds-1.);

  //   vector<double> bounds;
  //   bounds.push_back(atm.rmin);
  //   for (int i=1;i<n_bounds-1;i++) {
  //     double tau_target = taueff_max - i*taueff_step;
  //     auto tau_n = [=](double r) { return std::make_pair(taueff_spline(r) - tau_target, taueff_spline.prime(r)); };
  //     std::uintmax_t iterations = 1000;
  //     double rtarget = boost::math::tools::newton_raphson_iterate(tau_n,
  // 								  1.01*bounds[i-1], atm.rmin, atm.rmax,
  // 								  20, iterations);
  //     bounds.push_back(rtarget);
  //   }
  //   bounds.push_back(atm.rmax);

  //   return bounds;
  // }



//};








#endif
