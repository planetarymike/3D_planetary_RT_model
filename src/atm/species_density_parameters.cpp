#include "species_density_parameters.hpp"
using std::vector;
using std::exp;

#include <boost/numeric/odeint/integrate/integrate.hpp>
using boost::numeric::odeint::integrate;

#include <boost/numeric/odeint/integrate/integrate_const.hpp>
using boost::numeric::odeint::integrate_const;

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
using boost::numeric::odeint::runge_kutta4;

#include "push_back.hpp"

species_density_parameters::species_density_parameters(const double masss,
						       const double alphaa,
						       diffusion_coefs difff)
  : mass(masss), alpha(alphaa), diff(difff)
{ }

double species_density_parameters::get_CO2_exobase_density(const double &nCO2_lower,
							   const double &r_lower,
							   const double &rexoo,
							   temperature *tempp) {
  rexo = rexoo;
  temp = tempp;
  escape_flux = 0.0;
  
  double nCO2_exo, n_species_exo;
  get_thermosphere_density_exact(nCO2_lower,
				 /*n_species_lower*/ 1.0,
				 r_lower,
				 
				 escape_flux,
				 tempp,
				 
				 rexo,
				 nCO2_exo,
				 n_species_exo);

  return nCO2_exo;
}

// routine to get exact densities
void species_density_parameters::get_thermosphere_density_exact(const double &n_species_exo,
								const double &n_CO2_exo,
								const double &rexoo,
								 
								const double &escape_fluxx,
								temperature *tempp,
								 
								const double &r, 
								double &nCO2,
								double &n_species) {

  vector<double> log_nCO2_thermosphere_tmp;
  vector<double> log_n_species_thermosphere_tmp;
  vector<double> r_thermosphere_tmp;

  get_thermosphere_density_arrays(n_species_exo,
				  n_CO2_exo,
				  rexoo,
				  r,
				  
				  escape_fluxx,
				  tempp,
				  
				  log_nCO2_thermosphere_tmp,
				  log_n_species_thermosphere_tmp,
				  r_thermosphere_tmp,
				  200,
				  /*get_interpolation_points = */false);

  nCO2 = exp(log_nCO2_thermosphere_tmp.back());
  n_species = exp(log_n_species_thermosphere_tmp.back());
}



hydrogen_density_parameters::hydrogen_density_parameters()
  : species_density_parameters(mH, alpha_hydrogen, diffusion_coefs(DH0_hydrogen, s_hydrogen) )
{ }

void hydrogen_density_parameters::operator()( const vector<double> &x , vector<double> &dxdr , const double &r ) {
  // this operator returns the derivatives of log(nCO2) and log(nH)
  // for use with the Boost differential equations library
  
  //  x[0] = log(nCO2)
  //  x[1] = log(nH)
  
  // get the temperaure at this location
  const double temp_T      = temp->T(r);
  const double temp_Tprime = temp->Tprime(r);

  // set the diffusion coefficients (stored inside diff)
  diff.get(temp_T, temp->T_exo, exp(x[0]) );

  // get the inverse neutral (CO2) scale height
  double Hninv = G*mMars*mCO2/(kB*temp_T*r*r)+temp_Tprime/temp_T;

  // now get the inverse H scale height
  double HHinv = G*mMars*mass/(kB*temp_T*r*r)+(1+alpha)*temp_Tprime/temp_T;
  
  // now return the derivative of d(log(nCO2))/dr = - 1 / neutral scale height
  dxdr[0] = -Hninv; 

  // and d(log(nH))/dr = - 1/(D + K)*[ phi_H*(r_exo/r)^2/nH
  //                                   + D*(1 / H scale height) + K*(1 / neutral scale height) ]
  dxdr[1] = -1./(diff.DH+diff.KK)*( escape_flux * ( rexo*rexo /r /r ) / exp(x[1])
				    + (diff.DH*HHinv + diff.KK*Hninv) );
}

// called by thermosphere_exosphere to populate thermosphere interpolation arrays
void hydrogen_density_parameters::get_thermosphere_density_arrays(const double &n_species_exo,
								  const double &n_CO2_exo,
								  const double &rexoo,
								  const double &r, // minimum altitude to integrate to
								  
								  const double &escape_fluxx,
								  temperature *tempp,
								  
								  vector<double> &log_nCO2_thermosphere,
								  vector<double> &log_n_species_thermosphere,
								  vector<double> &r_thermosphere,
								  const int n_thermosphere_steps,
								  bool get_interpolation_points) {
  rexo = rexoo;
  escape_flux = escape_fluxx;
  temp = tempp;
  
  // reset the thermosphere vectors
  log_nCO2_thermosphere.clear();
  log_n_species_thermosphere.clear();
  r_thermosphere.clear();
  
  //integrate to get the species densities in the thermosphere
  vector<double> nexo(2);
  nexo[0] = log(n_CO2_exo);
  nexo[1] = log(n_species_exo);

  double thermosphere_step_r = -(rexo-r)/(n_thermosphere_steps-1.);

  if (get_interpolation_points) {
    //use a constant stepper for easy interpolation
    runge_kutta4< vector<double> > stepper;
    integrate_const( stepper , *this ,
		     nexo , rexo , r , thermosphere_step_r,
		     push_back_quantities( &log_nCO2_thermosphere,
					   &log_n_species_thermosphere,
					   &r_thermosphere ) );
  } else {
    // integrate with a variable step size to the requested r
    integrate( *this , nexo , rexo , r , thermosphere_step_r);
    log_nCO2_thermosphere.push_back(nexo[0]);
    log_n_species_thermosphere.push_back(nexo[1]);
  }
}


// atomic oxygen

oxygen_density_parameters::oxygen_density_parameters()
  : species_density_parameters(16*mH, alpha_oxygen, diffusion_coefs(DH0_oxygen, s_oxygen) )
{ }

void oxygen_density_parameters::operator()( const vector<double> &x , vector<double> &dxdr , const double &r ) {
  // this operator returns the derivatives of log(nCO2) and log(nH)
  // for use with the Boost differential equations library
  
  //  x[0] = log(nCO2)
  //  x[1] = log(nH)
  
  // get the temperaure at this location
  const double temp_T      = temp->T(r);
  const double temp_Tprime = temp->Tprime(r);

  // set the diffusion coefficients (stored inside diff)
  diff.get(temp_T, temp->T_exo, exp(x[0]) );

  // get the inverse neutral (CO2) scale height
  double Hninv = G*mMars*mCO2/(kB*temp_T*r*r)+temp_Tprime/temp_T;

  // now get the inverse H scale height
  double HHinv = G*mMars*mass/(kB*temp_T*r*r)+(1+alpha)*temp_Tprime/temp_T;
  
  // now return the derivative of d(log(nCO2))/dr = - 1 / neutral scale height
  dxdr[0] = -Hninv; 

  // co2 column
  double nCO2 = exp(x[0]);
  double overhead_co2 = exp(nCO2)/Hninv; // = n*H

  // determine the downward flux from the CO2 dissociation above this point
  double co2_N_ref = 2e16;
  double co2_b = -0.527553;
  double co2_jinf = 3.09e-7; // dissociation rate at 300 km
  // Justin determined these parameters, they are the best fit for the
  // expression j = J/nCO2 = j0*(1.0+N_co2/const)^power, which he fit
  // against my/Eryn's detailed photochemical model output for O.

  // We integrate the total dissociation rate of CO2 to get the
  // downward flux,

  // int_z0^inf J dz = int_z0^inf (J/nCO2) * (nCO2*dz)
  //
  //                 = int_0^N0 j dN,
  //
  // an integral we can perform analytically with respect to the total
  // column at this location:
  double flux_O_dissoc = co2_jinf * ((1.0 - std::pow(1.0 + std::abs(overhead_co2)/co2_N_ref, co2_b+1))
				     /
				     ((co2_b+1) / co2_N_ref));


  // to handle the situation in the lower atmosphere we must determine
  // the flux near photochemical equilibrium
  double O_recomb_rate = 1.8 * 3e-33 * std::pow(300.0 / temp_T, 3.25);
  double local_co2_dissoc_rate = co2_jinf*std::pow(1.0 + std::abs(overhead_co2)/co2_N_ref, co2_b);
  double O_equilibrium_density = std::sqrt(local_co2_dissoc_rate / (2.0*O_recomb_rate));
  double inv_J_scale_height = -co2_b*Hninv; // 1 / (photolysis scale height)

  double flux_O_near_equilibrium = -(diff.KK
				     *O_equilibrium_density
				     *(Hninv+0.5*inv_J_scale_height)
				     *(1+(1.25*Hninv
					  *0.5*(Hninv + inv_J_scale_height)
					  *(diff.KK*O_equilibrium_density)
					  /(local_co2_dissoc_rate*exp(x[0])))));
  // -(diff.KK * O_equilibrium_density * (Hninv + 0.5*inv_J_scale_height)
  // 				     / (1 + (1.25 * Hninv * 0.5*(Hninv + inv_J_scale_height)
  // 					     * (diff.KK * O_equilibrium_density
  // 						/ (local_co2_dissoc_rate * exp(x[0]))))));
  // this expression was derived by Justin by linearizing the
  // expression for vertical flux assuming a small departure from
  // photochemical equilibrium
  
  // to bridge the lower and upper atmosphere use the harmonic mean of both fluxes:
  // TODO: figure out what the problem is here
  double flux_O = 0.0;//1/(1/flux_O_dissoc + 1/flux_O_near_equilibrium);
  
  // and d(log(nO))/dr = - 1/(D + K)*[ phi_O*(r_exo/r)^2/nO
  //                                   - O_flux/nO
  //                                   + D*(1 / H scale height) + K*(1 / neutral scale height) ]
  double nO = exp(x[1]);
  dxdr[1] = -1./(diff.DH+diff.KK)*( escape_flux * ( rexo*rexo /r /r ) / nO
				    + flux_O / nO
				    + (diff.DH*HHinv + diff.KK*Hninv) );
  assert(!std::isnan(dxdr[1]) && "derivative values must be real numbers");
}

// called by thermosphere_exosphere to populate thermosphere interpolation arrays
void oxygen_density_parameters::get_thermosphere_density_arrays(const double &n_species_exo,
								const double &n_CO2_exo,
								const double &rexoo,
								const double &r, // minimum altitude to integrate to
								
								const double &escape_fluxx,
								temperature *tempp,
								
								vector<double> &log_nCO2_thermosphere,
								vector<double> &log_n_species_thermosphere,
								vector<double> &r_thermosphere,
								const int n_thermosphere_steps,
								bool get_interpolation_points) {
  rexo = rexoo;
  escape_flux = escape_fluxx;
  temp = tempp;
  
  // reset the thermosphere vectors
  log_nCO2_thermosphere.clear();
  log_n_species_thermosphere.clear();
  r_thermosphere.clear();
  
  //integrate to get the species densities in the thermosphere
  vector<double> nexo(2);
  nexo[0] = log(n_CO2_exo);
  nexo[1] = log(n_species_exo); // atomic oxygen density

  double thermosphere_step_r = -(rexo-r)/(n_thermosphere_steps-1.);

  if (get_interpolation_points) {
    //use a constant stepper for easy interpolation
    runge_kutta4< vector<double> > stepper;
    integrate_const( stepper , *this ,
		     nexo , rexo , r , thermosphere_step_r,
		     push_back_quantities( &log_nCO2_thermosphere,
					   &log_n_species_thermosphere,
					   &r_thermosphere ) );
  } else {
    // integrate with a variable step size to the requested r
    integrate( *this , nexo , rexo , r , thermosphere_step_r);
    log_nCO2_thermosphere.push_back(nexo[0]);
    log_n_species_thermosphere.push_back(nexo[2]);
  }
}
