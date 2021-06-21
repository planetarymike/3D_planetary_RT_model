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
						       const int n_varss,
						       diffusion_coefs difff)
  : mass(masss), alpha(alphaa), n_vars(n_varss), diff(difff)
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
  : species_density_parameters(mH, alpha_hydrogen, 2, diffusion_coefs(DH0_hydrogen, s_hydrogen) )
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
  vector<double> nexo(n_vars);
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
