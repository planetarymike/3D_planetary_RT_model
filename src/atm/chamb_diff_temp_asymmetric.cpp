#include "chamb_diff_temp_asymmetric.hpp"

doubReal chamb_diff_temp_asymmetric::T_sza(const doubReal &sza) const {
  return T0 + (T1-T0)*sza/pi;
}

doubReal chamb_diff_temp_asymmetric::n_species_sza(const doubReal sza) const {
  return n_species_sza(sza, A);
}

doubReal chamb_diff_temp_asymmetric::n_species_sza(const doubReal &sza, const doubReal AA/*=A*/) const {
  return AA*std::pow(T_sza(sza),-Tpower);
}

doubReal chamb_diff_temp_asymmetric::sza_int(const doubReal &f0, const doubReal &f1,
					   const doubReal &t0, const doubReal &t1) const {
  doubReal m = (f1-f0)/(t1-t0);

  using std::sin;
  using std::cos;

  return m*(sin(t1)-sin(t0))-(f1*cos(t1)-f0*cos(t0));
}

chamb_diff_temp_asymmetric::chamb_diff_temp_asymmetric(species_density_parameters *species_thermospheree,
						       const doubReal navgg,
						       const doubReal T00,
						       const doubReal T11)
  : chamb_diff_temp_asymmetric(species_thermospheree,
			       navgg,
			       T00,
			       T11,
			       /*      nCO2rmin = */2.6e13,			       
			       /*          rexo = */rexo_typical,
			       /*          rmin = */rMars + 80e5,
			       /*         rmaxx = */rMars + 50000e5,
			       /* rmindiffusion = */rMars + 80e5)
{ }

chamb_diff_temp_asymmetric::chamb_diff_temp_asymmetric(species_density_parameters *species_thermospheree,
						       const doubReal navgg,
						       const doubReal T00,
						       const doubReal T11,
						       const doubReal nCO2rminn, //a good number is 2.6e13 (80km) [Chaufray+2008]
						       const doubReal rexoo,
						       const doubReal rminn,
						       const doubReal rmaxx,
						       const doubReal rmindiffusionn,
						       //extra args for krasnopolsky_temp
						       const doubReal T_tropo/* = 125.0*/,
						       const doubReal r_tropo/* = rMars + 90e5*/,
						       const doubReal shape_parameter/* = 11.4*/,
						       //power for temperature in the expression n*T^p = const.
						       const doubReal Tpowerr/* = 2.5*/)
  : atmosphere(rminn, rexoo, rmaxx),
    navg(navgg), T0(T00), T1(T11), Tpower(Tpowerr), nCO2rmin(nCO2rminn), rmindiffusion(rmindiffusionn)
{
  //compute the normalizing constant for this T0 and T1
  const int n_sza_int=100;
  const doubReal dtheta = M_PI/n_sza_int;
  doubReal theta = 0.0;
  doubReal Atemp = 0.5*n_species_sza(theta, 1.0)*std::sin(theta); // I know it's zero but for completeness
  doubReal Anorm = 0.5*std::sin(theta);
  for (int i=1;i<n_sza_int;i++) {
    theta+=dtheta;
    Atemp+=n_species_sza(theta,1.0)*std::sin(theta);
    Anorm+=std::sin(theta);
  }
  theta+=dtheta;//now theta = pi
  Atemp+=0.5*n_species_sza(theta,1.0)*std::sin(theta);
  Anorm+=0.5*std::sin(theta);
  
  A = navg*Anorm/Atemp;//now we can call n_species_sza(sza) and get the correct result

  //set up each of the temperatures and atmospheres
  sza_vec.resize(n_sza);
  for (int isza=0; isza<n_sza; isza++) {
    sza_vec[isza] = isza*d_sza;
    Temp_sza[isza] = krasnopolsky_temperature(T_sza(sza_vec[isza]),
					      T_tropo,
					      r_tropo,
					      shape_parameter,
					      false/*shape_parameter is in absolute units of km*/);

    atm_sza[isza] = new chamb_diff_1d(rmin,
				      rexo,
				      rmax,
				      rmindiffusion,
				      n_species_sza(sza_vec[isza]),
				      nCO2rmin,
				      &Temp_sza[isza],
				      species_thermospheree,
				      thermosphere_exosphere::method_rmax_nCO2rmin);
  }

  // //set the max altitude to the minumum of the max altitudes as a function of SZA
  // rmax = atm_sza[n_sza-1]->rmax;
  // for (int isza=0; isza<n_sza; isza++) {
  //   atm_sza[isza]->rmax = rmax;
  //   atm_sza[isza]->atmosphere_average_1d::setup();//redo vertical integration with new rmax
  // }
  //all of the 1d atmospheres use the same radial grid for averaging
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

  n_species_int_interp = Bilinear_interp<doubReal>(log_r,sza_vec,n_species_int);
  n_absorber_int_interp = Bilinear_interp<doubReal>(log_r,sza_vec,n_absorber_int);
  Temp_int_interp = Bilinear_interp<doubReal>(log_r,sza_vec,Temp_int);

  init=true;
}
  
chamb_diff_temp_asymmetric::~chamb_diff_temp_asymmetric() {
  //free allocated atmosphere memory
  for (int isza=0; isza<n_sza; isza++)
    delete atm_sza[isza];
}


void chamb_diff_temp_asymmetric::sza_interp(const doubReal &sza, int &i_sza, doubReal &sza_wt) const {
  i_sza = (int) (sza/d_sza);//rounds down (truncation)
  sza_wt = 1.0-(sza-sza_vec[i_sza])/d_sza;// 1.0 if we're at the lower boundary, 0.0 at the upper
}

doubReal chamb_diff_temp_asymmetric::r_to_log_r(const doubReal &r) const {
  return std::log((r-rMars)/atmosphere_average_1d::r_int_scale);
}

doubReal chamb_diff_temp_asymmetric::avg(const Bilinear_interp<doubReal> &terp,
				       const doubReal &r0, const doubReal &r1,
				       const doubReal &t0, const doubReal &t1) const {
  //     q0   |    |
  //          | q1 |
  //  --------x----x
  //     q2   | q3 |
  //  --------x----x
  
  doubReal q0 = terp.interp(r_to_log_r(r1),t0);
  doubReal q1 = terp.interp(r_to_log_r(r1),t1);
  doubReal q2 = terp.interp(r_to_log_r(r0),t0);
  doubReal q3 = terp.interp(r_to_log_r(r0),t1);

  doubReal qint = q3 - q2 - q1 + q0;

  //need to use the same scaling as atmosphere_average_1d to keep the units correct
  doubReal r00 = r0/atmosphere_average_1d::r_int_scale;
  doubReal r11 = r1/atmosphere_average_1d::r_int_scale;

  doubReal rweight = (r11*r11*r11-r00*r00*r00)/3.0;// int_r0^r1 (r^2) dr
  doubReal szaweight = std::cos(t0) - std::cos(t1);// int_t0^t1 sin(t) dt

  return qint/rweight/szaweight;
}



// doubReal chamb_diff_temp_asymmetric::H_Temp(const atmo_point &pt) const {
//   return Temp(pt);
// }
// void chamb_diff_temp_asymmetric::H_Temp(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
//   Temp(vox, ret_avg, ret_pt);
// }
doubReal chamb_diff_temp_asymmetric::Temp(const atmo_point &pt) const {
  if (!temp_dependent_sH)
    return constant_temp_sH;
  else {
    int isza;
    doubReal szawt;
    sza_interp(pt.t, isza, szawt);
    if (isza==n_sza-1)
      return atm_sza[isza]->thermosphere_exosphere::Temp(pt.r);
    else
      return (       szawt *(atm_sza[isza  ]->thermosphere_exosphere::Temp(pt.r))
	      + (1.0-szawt)*(atm_sza[isza+1]->thermosphere_exosphere::Temp(pt.r)));
  }
}
void chamb_diff_temp_asymmetric::Temp_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  if (!temp_dependent_sH)
    ret_avg = ret_pt = constant_temp_sH;
  else {
    doubReal szamin = vox.tbounds[0];
    if (szamin < 0)
      szamin=0;
    doubReal szamax = vox.tbounds[1];
    if (szamax>M_PI)
      szamax = M_PI;

    // //test points at the corners to validate averages
    // doubReal test_avg[4];
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

doubReal chamb_diff_temp_asymmetric::n_absorber(const atmo_point &pt) const {
  return nCO2(pt);
} 
void chamb_diff_temp_asymmetric::n_absorber_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  nCO2_voxel_avg(vox, ret_avg, ret_pt);
}
doubReal chamb_diff_temp_asymmetric::nCO2(const atmo_point &pt) const {
  int isza;
  doubReal szawt;
  sza_interp(pt.t, isza, szawt);
  if (isza==n_sza-1)
    return atm_sza[isza]->thermosphere_exosphere::nCO2(pt.r);
  else
    return (       szawt *(atm_sza[isza  ]->thermosphere_exosphere::nCO2(pt.r))
	    + (1.0-szawt)*(atm_sza[isza+1]->thermosphere_exosphere::nCO2(pt.r)));
}
void chamb_diff_temp_asymmetric::nCO2_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  doubReal szamin = vox.tbounds[0];
  if (szamin < 0)
    szamin=0;
  doubReal szamax = vox.tbounds[1];
  if (szamax>M_PI)
    szamax = M_PI;

  // //test points at the corners to validate averages
  // doubReal test_avg[4];
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

doubReal chamb_diff_temp_asymmetric::n_species(const atmo_point &pt) const {
  int isza;
  doubReal szawt;
  sza_interp(pt.t, isza, szawt);
  if (isza==n_sza-1)
    return atm_sza[isza]->thermosphere_exosphere::n_species(pt.r);
  else
    return (       szawt *(atm_sza[isza  ]->thermosphere_exosphere::n_species(pt.r))
	    + (1.0-szawt)*(atm_sza[isza+1]->thermosphere_exosphere::n_species(pt.r)));
}
void chamb_diff_temp_asymmetric::n_species_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  doubReal szamin = vox.tbounds[0];
  if (szamin < 0)
    szamin=0;
  doubReal szamax = vox.tbounds[1];
  if (szamax>M_PI)
    szamax = M_PI;

  // //test points at the corners to validate averages
  // doubReal test_avg[4];
  // atmo_point my_pt;
  // my_pt.rtp(vox.rbounds[0],szamin,0);
  // test_avg[0] = n_species(my_pt);
  // my_pt.rtp(vox.rbounds[1],szamin,0);
  // test_avg[1] = n_species(my_pt);
  // my_pt.rtp(vox.rbounds[0],szamax,0);
  // test_avg[2] = n_species(my_pt);
  // my_pt.rtp(vox.rbounds[1],szamax,0);
  // test_avg[3] = n_species(my_pt);

  ret_avg = avg( n_species_int_interp,
		 vox.rbounds[0], vox.rbounds[1],
		 szamin        , szamax);

  ret_pt = n_species(vox.pt);
}

doubReal chamb_diff_temp_asymmetric::n_species(const doubReal &r) const {
  return atm_sza[0]->thermosphere_exosphere::n_species(r);
}

doubReal chamb_diff_temp_asymmetric::r_from_n_species(const doubReal &n_species) const {
  return atm_sza[0]->r_from_n_species(n_species);
}

doubReal chamb_diff_temp_asymmetric::Temp(const doubReal &r) const {
  return atm_sza[0]->thermosphere_exosphere::Temp(r);
}

doubReal chamb_diff_temp_asymmetric::n_absorber(const doubReal &r) const {
  return atm_sza[0]->thermosphere_exosphere::n_absorber(r);
}
