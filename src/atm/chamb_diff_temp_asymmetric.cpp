#include "chamb_diff_temp_asymmetric.hpp"

double chamb_diff_temp_asymmetric::T_sza(const double &sza) const {
  return T0 + (T1-T0)*sza/pi;
}

double chamb_diff_temp_asymmetric::nH_sza(const double sza) const {
  return nH_sza(sza,A);
}

double chamb_diff_temp_asymmetric::nH_sza(const double &sza, const double AA/*=A*/) const {
  return AA*std::pow(T_sza(sza),-Tpower);
}

double chamb_diff_temp_asymmetric::sza_int(const double &f0, const double &f1,
					   const double &t0, const double &t1) const {
  double m = (f1-f0)/(t1-t0);

  using std::sin;
  using std::cos;

  return m*(sin(t1)-sin(t0))-(f1*cos(t1)-f0*cos(t0));
}

chamb_diff_temp_asymmetric::chamb_diff_temp_asymmetric(const double navgg,
						       const double T00,
						       const double T11)
  : chamb_diff_temp_asymmetric(navgg,
			       T00,
			       T11,
			       /*      nCO2rmin = */2.6e13,			       
			       /*          rexo = */rexo_typical,
			       /*          rmin = */rMars + 80e5,
			       /*         rmaxx = */rMars + 50000e5,
			       /* rmindiffusion = */rMars + 80e5)
{ }

chamb_diff_temp_asymmetric::chamb_diff_temp_asymmetric(const double navgg,
						       const double T00,
						       const double T11,
						       const double nCO2rminn, //a good number is 2.6e13 (80km) [Chaufray+2008]
						       const double rexoo,
						       const double rminn,
						       const double rmaxx,
						       const double rmindiffusionn,
						       //extra args for krasnopolsky_temp
						       const double T_tropo/* = 125.0*/,
						       const double r_tropo/* = rMars + 90e5*/,
						       const double shape_parameter/* = 11.4*/)
  : atmosphere(rminn, rexoo, rmaxx),
    navg(navgg), T0(T00), T1(T11), nCO2rmin(nCO2rminn), rmindiffusion(rmindiffusionn)
{
  //compute the normalizing constant for this T0 and T1
  const int n_sza_int=100;
  const double dtheta = M_PI/n_sza_int;
  double theta = 0.0;
  double Atemp = 0.5*nH_sza(theta, 1.0)*std::sin(theta); // I know it's zero but for completeness
  double Anorm = 0.5*std::sin(theta);
  for (int i=1;i<n_sza_int;i++) {
    theta+=dtheta;
    Atemp+=nH_sza(theta,1.0)*std::sin(theta);
    Anorm+=std::sin(theta);
  }
  theta+=dtheta;//now theta = pi
  Atemp+=0.5*nH_sza(theta,1.0)*std::sin(theta);
  Anorm+=0.5*std::sin(theta);
  
  A = navg*Anorm/Atemp;//now we can call nH_sza(sza) and get the correct result

  //set up each of the temperatures and atmospheres
  sza_vec.resize(n_sza);
  for (int isza=0; isza<n_sza; isza++) {
    sza_vec[isza] = isza*d_sza;
    Temp_sza[isza] = krasnopolsky_temperature(T_sza(sza_vec[isza]),
					      T_tropo,
					      r_tropo,
					      shape_parameter);

    atm_sza[isza] = new chamb_diff_1d(rmin,
				      rexo,
				      rmax,
				      rmindiffusion,
				      nH_sza(sza_vec[isza]),
				      nCO2rmin,
				      Temp_sza[isza],
				      thermosphere_exosphere::method_rmax_nCO2rmin);
  }

  //set the max altitude to the minumum of the max altitudes as a function of SZA
  rmax = atm_sza[n_sza-1]->rmax;
  for (int isza=0; isza<n_sza; isza++) {
    atm_sza[isza]->rmax = rmax;
    atm_sza[isza]->atmosphere_average_1d::setup();//redo vertical integration with new rmax
  }
  //now all of the 1d atmospheres use the same radial grid for averaging
  log_r = atm_sza[n_sza-1]->log_r_int;

  
  //each of the 1d atmospheres has already been integrated in radius
  //now we need to integrate over sza to get 2d averages
  n_species_int.resize(n_log_r, n_sza);
  n_absorber_int.resize(n_log_r, n_sza);
  Temp_int.resize(n_log_r, n_sza);
  for (int ir = 0; ir<n_log_r; ir++) {
    n_species_int(ir,0) = 0.0;
    n_absorber_int(ir,0) = 0.0;
    Temp_int(ir,0) = 0.0;
    for (int isza = 1; isza < n_sza; isza++) {
      n_species_int(ir,isza)  = n_species_int(ir,isza-1) +
                                 sza_int(atm_sza[isza - 1]->n_species_int_spherical[ir], atm_sza[isza]->n_species_int_spherical[ir],
                                         sza_vec[isza - 1]                             , sza_vec[isza]);
      n_absorber_int(ir,isza) = n_absorber_int(ir,isza-1) +
                                 sza_int(atm_sza[isza - 1]->n_absorber_int_spherical[ir], atm_sza[isza]->n_absorber_int_spherical[ir],
                                         sza_vec[isza - 1]                              , sza_vec[isza]);
      Temp_int(ir,isza)       = Temp_int(ir,isza-1) +
	                         sza_int(atm_sza[isza - 1]->Tint_spherical[ir], atm_sza[isza]->Tint_spherical[ir],
					 sza_vec[isza - 1]                    , sza_vec[isza]);
    }
  }
  //integration is done. Now we can create 2D interpolation objects to read average densities from

  n_species_int_interp = Bilinear_interp<double>(log_r,sza_vec,n_species_int);
  n_absorber_int_interp = Bilinear_interp<double>(log_r,sza_vec,n_absorber_int);
  Temp_int_interp = Bilinear_interp<double>(log_r,sza_vec,Temp_int);

  init=true;
}
  
chamb_diff_temp_asymmetric::~chamb_diff_temp_asymmetric() {
  //free allocated atmosphere memory
  for (int isza=0; isza<n_sza; isza++)
    delete atm_sza[isza];
}


void chamb_diff_temp_asymmetric::sza_interp(const double &sza, int &i_sza, double &sza_wt) const {
  i_sza = (int) (sza/d_sza);//rounds down (truncation)
  sza_wt = 1.0-(sza-sza_vec[i_sza])/d_sza;// 1.0 if we're at the lower boundary, 0.0 at the upper
}

double chamb_diff_temp_asymmetric::r_to_log_r(const double &r) const {
  return std::log((r-rMars)/atmosphere_average_1d::r_int_scale);
}

double chamb_diff_temp_asymmetric::avg(const Bilinear_interp<double> &terp,
				       const double &r0, const double &r1,
				       const double &t0, const double &t1) const {
  //     q0   |    |
  //          | q1 |
  //  --------x----x
  //     q2   | q3 |
  //  --------x----x
  
  double q0 = terp.interp(r_to_log_r(r1),t0);
  double q1 = terp.interp(r_to_log_r(r1),t1);
  double q2 = terp.interp(r_to_log_r(r0),t0);
  double q3 = terp.interp(r_to_log_r(r0),t1);

  double qint = q3 - q2 - q1 + q0;

  //need to use the same scaling as atmosphere_average_1d to keep the units correct
  double r00 = r0/atmosphere_average_1d::r_int_scale;
  double r11 = r1/atmosphere_average_1d::r_int_scale;

  double rweight = (r11*r11*r11-r00*r00*r00)/3.0;// int_r0^r1 (r^2) dr
  double szaweight = std::cos(t0) - std::cos(t1);// int_t0^t1 sin(t) dt

  return qint/rweight/szaweight;
}



Real chamb_diff_temp_asymmetric::H_Temp(const atmo_point &pt) const {
  return Temp(pt);
}
void chamb_diff_temp_asymmetric::H_Temp(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  Temp(vox, ret_avg, ret_pt);
}
Real chamb_diff_temp_asymmetric::Temp(const atmo_point &pt) const {
  if (!temp_dependent_sH)
    return constant_temp_sH;
  else {
    int isza;
    double szawt;
    sza_interp(pt.t, isza, szawt);
    if (isza==n_sza-1)
      return atm_sza[isza]->thermosphere_exosphere::Temp(pt.r);
    else
      return (       szawt *(atm_sza[isza  ]->thermosphere_exosphere::Temp(pt.r))
	      + (1.0-szawt)*(atm_sza[isza+1]->thermosphere_exosphere::Temp(pt.r)));
  }
}
void chamb_diff_temp_asymmetric::Temp(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  if (!temp_dependent_sH)
    ret_avg = ret_pt = constant_temp_sH;
  else {
    double szamin = vox.tbounds[0];
    if (szamin < 0)
      szamin=0;
    double szamax = vox.tbounds[1];
    if (szamax>M_PI)
      szamax = M_PI;

    // //test points at the corners to validate averages
    // double test_avg[4];
    // atmo_point my_pt;
    // my_pt.rtp(vox.rbounds[0],szamin,0);
    // test_avg[0] = Temp(my_pt);
    // my_pt.rtp(vox.rbounds[1],szamin,0);
    // test_avg[1] = Temp(my_pt);
    // my_pt.rtp(vox.rbounds[0],szamax,0);
    // test_avg[2] = Temp(my_pt);
    // my_pt.rtp(vox.rbounds[1],szamax,0);
    // test_avg[3] = Temp(my_pt);

    ret_avg = avg( Temp_int_interp,
		   vox.rbounds[0], vox.rbounds[1],
		   szamin        , szamax);
    
    ret_pt = Temp(vox.pt);
  }
}

Real chamb_diff_temp_asymmetric::n_absorber(const atmo_point &pt) const {
  return nCO2(pt);
} 
void chamb_diff_temp_asymmetric::n_absorber(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  nCO2(vox, ret_avg, ret_pt);
}
Real chamb_diff_temp_asymmetric::nCO2(const atmo_point &pt) const {
  int isza;
  double szawt;
  sza_interp(pt.t, isza, szawt);
  if (isza==n_sza-1)
    return atm_sza[isza]->thermosphere_exosphere::nCO2(pt.r);
  else
    return (       szawt *(atm_sza[isza  ]->thermosphere_exosphere::nCO2(pt.r))
	    + (1.0-szawt)*(atm_sza[isza+1]->thermosphere_exosphere::nCO2(pt.r)));
}
void chamb_diff_temp_asymmetric::nCO2(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  double szamin = vox.tbounds[0];
  if (szamin < 0)
    szamin=0;
  double szamax = vox.tbounds[1];
  if (szamax>M_PI)
    szamax = M_PI;

  // //test points at the corners to validate averages
  // double test_avg[4];
  // atmo_point my_pt;
  // my_pt.rtp(vox.rbounds[0],szamin,0);
  // test_avg[0] = nCO2(my_pt);
  // my_pt.rtp(vox.rbounds[1],szamin,0);
  // test_avg[1] = nCO2(my_pt);
  // my_pt.rtp(vox.rbounds[0],szamax,0);
  // test_avg[2] = nCO2(my_pt);
  // my_pt.rtp(vox.rbounds[1],szamax,0);
  // test_avg[3] = nCO2(my_pt);
  
  ret_avg = avg( n_absorber_int_interp,
		 vox.rbounds[0], vox.rbounds[1],
		 szamin        , szamax);

  ret_pt = nCO2(vox.pt);
}

Real chamb_diff_temp_asymmetric::n_species(const atmo_point &pt) const {
  return nH(pt);
}
void chamb_diff_temp_asymmetric::n_species(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  nH(vox, ret_avg, ret_pt);
}
Real chamb_diff_temp_asymmetric::nH(const atmo_point &pt) const {
  int isza;
  double szawt;
  sza_interp(pt.t, isza, szawt);
  if (isza==n_sza-1)
    return atm_sza[isza]->thermosphere_exosphere::nH(pt.r);
  else
    return (       szawt *(atm_sza[isza  ]->thermosphere_exosphere::nH(pt.r))
	    + (1.0-szawt)*(atm_sza[isza+1]->thermosphere_exosphere::nH(pt.r)));
}
void chamb_diff_temp_asymmetric::nH(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  double szamin = vox.tbounds[0];
  if (szamin < 0)
    szamin=0;
  double szamax = vox.tbounds[1];
  if (szamax>M_PI)
    szamax = M_PI;

  // //test points at the corners to validate averages
  // double test_avg[4];
  // atmo_point my_pt;
  // my_pt.rtp(vox.rbounds[0],szamin,0);
  // test_avg[0] = nH(my_pt);
  // my_pt.rtp(vox.rbounds[1],szamin,0);
  // test_avg[1] = nH(my_pt);
  // my_pt.rtp(vox.rbounds[0],szamax,0);
  // test_avg[2] = nH(my_pt);
  // my_pt.rtp(vox.rbounds[1],szamax,0);
  // test_avg[3] = nH(my_pt);

  ret_avg = avg( n_species_int_interp,
		 vox.rbounds[0], vox.rbounds[1],
		 szamin        , szamax);

  ret_pt = nH(vox.pt);
}

Real chamb_diff_temp_asymmetric::n_species(const Real &r) const {
  return atm_sza[0]->thermosphere_exosphere::n_species(r);
}

Real chamb_diff_temp_asymmetric::r_from_n_species(const Real &n_species) const {
  return atm_sza[0]->r_from_n_species(n_species);
}

Real chamb_diff_temp_asymmetric::Temp(const Real &r) const {
  return atm_sza[0]->thermosphere_exosphere::Temp(r);
}

Real chamb_diff_temp_asymmetric::n_absorber(const Real &r) const {
  return atm_sza[0]->thermosphere_exosphere::n_absorber(r);
}
