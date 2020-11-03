//lineshape_tracker.cpp --- holstein integrals computed JIT as lines of sight are traversed
#include "constants.hpp"
#include "lineshape_tracker.hpp"

// CUDA_CALLABLE_MEMBER
// lineshape_tracker::lineshape_tracker() { }
// CUDA_CALLABLE_MEMBER
// lineshape_tracker::~lineshape_tracker() { }
// CUDA_CALLABLE_MEMBER
// lineshape_tracker::lineshape_tracker(const lineshape_tracker &copy)
// {
//   //lambda = {SHIZGAL_LAMBDAS};
//   //weight = {SHIZGAL_WEIGHTS};
  
//   for (int i_lambda = 0; i_lambda<N_LAMBDA; i_lambda++) {
//     // lambda[i_lambda]   = copy.lambda[i_lambda];
//     // weight[i_lambda]   = copy.weight[i_lambda];

//     //lambda2[i_lambda]  = copy.lambda2[i_lambda];
//     //weightfn[i_lambda] = copy.weightfn[i_lambda];

//     lineshape_at_origin[i_lambda] = copy.lineshape_at_origin[i_lambda];
//     //lineshape[i_lambda] = copy.lineshape[i_lambda];

//     tau_species_lambda_initial[i_lambda] = copy.tau_species_lambda_initial[i_lambda];
//     //tau_species_lambda_final[i_lambda] = copy.tau_species_lambda_final[i_lambda];

//     transfer_probability_lambda_initial[i_lambda] = copy.transfer_probability_lambda_initial[i_lambda];
//     // transfer_probability_lambda_voxel[i_lambda] = copy.transfer_probability_lambda_voxel[i_lambda];
//     // transfer_probability_lambda_final[i_lambda] = copy.transfer_probability_lambda_final[i_lambda];
//   }

//   tau_species_initial=copy.tau_species_initial;
//   tau_species_final=copy.tau_species_final;

//   tau_absorber_initial=copy.tau_absorber_initial;
//   tau_absorber_final=copy.tau_absorber_final;

//   max_tau_species=copy.max_tau_species;

//   holstein_T_initial = copy.holstein_T_final;
//   holstein_T_final = copy.holstein_T_final;
//   holstein_T_int = copy.holstein_T_int;
//   holstein_G_int = copy.holstein_G_int;
// }
// CUDA_CALLABLE_MEMBER
// lineshape_tracker& lineshape_tracker::operator=(const lineshape_tracker &rhs) {
//   if(this == &rhs) return *this;

//   for (int i_lambda = 0; i_lambda<N_LAMBDA; i_lambda++) {
//     //lambda2[i_lambda]  = rhs.lambda2[i_lambda];
//     //weightfn[i_lambda] = rhs.weightfn[i_lambda];

//     lineshape_at_origin[i_lambda] = rhs.lineshape_at_origin[i_lambda];
//     //lineshape[i_lambda] = rhs.lineshape[i_lambda];

//     tau_species_lambda_initial[i_lambda] = rhs.tau_species_lambda_initial[i_lambda];
//     //tau_species_lambda_final[i_lambda] = rhs.tau_species_lambda_final[i_lambda];

//     transfer_probability_lambda_initial[i_lambda] = rhs.transfer_probability_lambda_initial[i_lambda];
//     // transfer_probability_lambda_voxel[i_lambda] = rhs.transfer_probability_lambda_voxel[i_lambda];
//     // transfer_probability_lambda_final[i_lambda] = rhs.transfer_probability_lambda_final[i_lambda];
//   }

//   tau_species_initial=rhs.tau_species_initial;
//   tau_species_final=rhs.tau_species_final;

//   tau_absorber_initial=rhs.tau_absorber_initial;
//   tau_absorber_final=rhs.tau_absorber_final;

//   max_tau_species=rhs.max_tau_species;

//   holstein_T_initial = rhs.holstein_T_final;
//   holstein_T_final = rhs.holstein_T_final;
//   holstein_T_int = rhs.holstein_T_int;
//   holstein_G_int = rhs.holstein_G_int;

//   return *this;
// }

CUDA_CALLABLE_MEMBER
Real lineshape_tracker::lambda(const int &i_lambda) const {
  return i_lambda*LAMBDA_MAX/(N_LAMBDA-1);
}
CUDA_CALLABLE_MEMBER
Real lineshape_tracker::weight(const int &i_lambda) const {
  Real thisweight = LAMBDA_MAX/(N_LAMBDA-1);
  if (i_lambda==0 || i_lambda==N_LAMBDA-1)
    thisweight *= REAL(0.5);
  return thisweight;
}
CUDA_CALLABLE_MEMBER
Real lineshape_tracker::line_shape_function(const Real &lambda2, const Real &T_ratio) const {
  return exp(-lambda2*T_ratio);
}
CUDA_CALLABLE_MEMBER
Real lineshape_tracker::line_shape_function_normalized(const Real &lambda2, const Real &T_ratio) const {
  return std::sqrt(T_ratio)*line_shape_function(lambda2,T_ratio);
}

CUDA_CALLABLE_MEMBER
void lineshape_tracker::init() {
  max_tau_species = 0;
}

CUDA_CALLABLE_MEMBER
void lineshape_tracker::reset(const Real &T_ratio) {
  tau_species_final=0;
  tau_absorber_final=0;
  holstein_T_final=1.0;
  T_ratio_at_origin=T_ratio;

  for (int i_lambda = 0; i_lambda<N_LAMBDA; i_lambda++) {
    //lineshape_at_origin[i_lambda] = std::sqrt(T_ratio)*exp(-lambda2*T_ratio);
    // assert(!std::isnan(lineshape_at_origin[i_lambda])
    // 	   && lineshape_at_origin[i_lambda]>0 && "lineshape must be real and positive");

    // tau_species_lambda_initial[i_lambda] = 0.0;
    transfer_probability_lambda_initial[i_lambda] = 1.0;
  }
  //max_tau_species not reset because we want to track this across
  //all lines of sight
}

CUDA_CALLABLE_MEMBER
void lineshape_tracker::check_max_tau() {
  //we only need to check the line center where the optical depth is greatest
  if (tau_species_final > max_tau_species)
    max_tau_species = tau_species_final;
}
  
CUDA_CALLABLE_MEMBER
void lineshape_tracker::update_start(const Real &T_ratio,
				     const Real &dtau_species,
				     const Real &dtau_absorber,
				     const Real &abs,
				     const Real &pathlength,
				     const bool influence/* = true*/)
{
  Real tau_species_voxel = dtau_species * pathlength;
  tau_species_final += tau_species_voxel;
  assert(!std::isnan(tau_species_final)
	 && tau_species_final>=0
	 && "optical depths must be real numbers");

  tau_absorber_final += dtau_absorber * pathlength; 
  assert(!std::isnan(tau_absorber_final)
	 && tau_absorber_final>=0
	 && "optical depths must be real numbers");
    

  holstein_T_final = 0;
  holstein_T_int = 0;
  holstein_G_int = 0;
  Real holTcoef;
  Real holstein_T_int_coef;

  for (int i_lambda = 0; i_lambda < N_LAMBDA; i_lambda++) {
    Real lambda2 = lambda(i_lambda) * lambda(i_lambda);
    Real lineshape = line_shape_function(lambda2, T_ratio);
    assert(!std::isnan(lineshape) && lineshape > 0 &&
           "lineshape must be real and positive");
    Real lineshape_at_origin = line_shape_function_normalized(lambda2, T_ratio_at_origin);
    assert(!std::isnan(lineshape_at_origin) && lineshape_at_origin>0 &&
	   "lineshape must be real and positive");

    // Real tau_species_lambda_final = (tau_species_lambda_initial[i_lambda] +
    //                                  (dtau_species * pathlength * lineshape));
    // assert(!std::isnan(tau_species_lambda_final) &&
    //        tau_species_lambda_final >= 0 &&
    //        "optical depth must be real and positive");

    Real transfer_probability_lambda_voxel = std::exp(-(dtau_absorber + dtau_species * lineshape) * pathlength);
    assert(!std::isnan(transfer_probability_lambda_voxel) &&
           transfer_probability_lambda_voxel >= 0 &&
           transfer_probability_lambda_voxel <= 1 &&
           "transfer probability is a probability.");

    Real transfer_probability_lambda_final = (transfer_probability_lambda_initial[i_lambda]
					      * transfer_probability_lambda_voxel);

    holTcoef = (two_over_sqrt_pi * weight(i_lambda));
    //		/ weightfn[i_lambda]);

    // holstein_T_int represents a frequency averaged emission and
    // therefore uses lineshape in this voxel
    holstein_T_int_coef = (holTcoef
			   * lineshape
			   * transfer_probability_lambda_initial[i_lambda]
			   * (REAL(1.0) - transfer_probability_lambda_voxel)
			   / (abs + lineshape));
    holstein_T_int += holstein_T_int_coef;
    assert(!std::isnan(holstein_T_int) && holstein_T_int >= 0 &&
           (holstein_T_int <= tau_species_voxel ||
            std::abs(holstein_T_int - tau_species_voxel) < ABS)
           //  ^^ this allows for small rounding errors
           && "holstein integral must be between 0 and Delta tau b/c 0<=HolT<=1");

    if (influence) {
      // holstein T final represents a frequency-averaged absorption and
      // uses lineshape_at_origin
      holstein_T_final += (holTcoef * lineshape_at_origin *
			   transfer_probability_lambda_final);
      assert(!std::isnan(holstein_T_final) && holstein_T_final >= 0 &&
	     holstein_T_final <= 1 &&
	     "holstein function represents a probability");

      // holstein_G_int represents the frequency averaged emission
      // followed by absorption and uses both lineshape in this voxel
      //(via holstein_T_int_coef) and the lineshape at origin
      holstein_G_int += holstein_T_int_coef * lineshape_at_origin;
      assert(!std::isnan(holstein_G_int) && holstein_G_int >= 0 && holstein_G_int <= 1 &&
	     "holstein G integral represents a probability");
    }

    //update the initial values
    //tau_species_lambda_initial[i_lambda] = tau_species_lambda_final;
    transfer_probability_lambda_initial[i_lambda] = transfer_probability_lambda_final;
  }
  //check that the last element is not contributing too much to the integral
  assert(!((holstein_T_int > STRICTABS) && (holstein_T_int_coef/holstein_T_int > 1e-2))
	 && "wings of line contribute too much to transmission. Increase LAMBDA_MAX in lineshape_tracker.");
  //if holstein T is larger than physically possible due to rounding errors, reduce it to the physical limit
  if (holstein_T_int > tau_species_voxel)
    holstein_T_int = tau_species_voxel; //we don't need to worry about
					//big rounding errors here
					//because we checked earlier
}
CUDA_CALLABLE_MEMBER
void lineshape_tracker::update_start_brightness(const Real &T_ratio,
						const Real &dtau_species,
						const Real &dtau_absorber,
						const Real &abs,
						const Real &pathlength)
{
  update_start(T_ratio,
	       dtau_species,
	       dtau_absorber,
	       abs,
	       pathlength,
	       /*influence = */false);
}

CUDA_CALLABLE_MEMBER
void lineshape_tracker::update_end() {
  // tau_species_initial=tau_species_final;
  // tau_absorber_initial=tau_absorber_final;
  // holstein_T_initial=holstein_T_final;
  
  // for (int i_lambda=0; i_lambda < N_LAMBDA; i_lambda++) {
  //   tau_species_lambda_initial[i_lambda] = tau_species_lambda_final[i_lambda];
  //   transfer_probability_lambda_initial[i_lambda] = transfer_probability_lambda_final[i_lambda];
  // }
  
  check_max_tau();
}

