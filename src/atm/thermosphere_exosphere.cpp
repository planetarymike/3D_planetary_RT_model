#include "thermosphere_exosphere.hpp"
using std::vector;

#include <fstream>

thermosphere_exosphere::thermosphere_exosphere(doubReal n_species_exoo, 
					       doubReal nCO2exoo, 
					       temperature *tempp,
					       species_density_parameters *species_thermospheree)
  : thermosphere_exosphere(/*          rmin = */rMars + 80e5,
			   /*          rexo = */rexo_typical,
			   /* n_species_min = */10,
			   /* rmindiffusion = */rMars + 80e5,
			   n_species_exoo,
			   nCO2exoo,
			   tempp,
			   species_thermospheree)   { }

thermosphere_exosphere::thermosphere_exosphere(doubReal rminn,
					       doubReal rexoo,
					       doubReal rmaxx_or_nspmin,
					       doubReal rmindiffusionn,
					       doubReal n_species_exoo, 
					       doubReal nCO2rmin_or_nCO2exoo, 
					       temperature *tempp,
					       species_density_parameters *species_thermospheree,
					       const int method)
  : atmosphere(rminn, rexoo, -1), // -1 is immediately overwritten
    exosphere(rexoo,
	      tempp->T_exo,
	      n_species_exoo,
	      species_thermospheree->mass), // overwritten in setup()
    species_thermosphere(species_thermospheree),
    CO2_exosphere(rexoo,
		  tempp->T_exo,
		  n_species_exoo,
		  mCO2)
{
  if (method == method_nspmin_nCO2exo)
    this->setup_nspmin_nCO2exo(rminn,
			       rexoo,
			       rmaxx_or_nspmin,
			       rmindiffusionn,
			       n_species_exoo,
			       nCO2rmin_or_nCO2exoo,
			       tempp);
  else if (method == method_rmax_nCO2rmin)
    this->setup_rmax_nCO2rmin(rminn,
			      rexoo,
			      rmaxx_or_nspmin,
			      rmindiffusionn,
			      n_species_exoo,
			      nCO2rmin_or_nCO2exoo,
			      tempp);
}

void thermosphere_exosphere::setup_nspmin_nCO2exo(doubReal rminn,
						 doubReal rexoo,
						 doubReal n_species_min,
						 doubReal rmindiffusionn,
						 doubReal n_species_exoo, 
						 doubReal nCO2exoo, 
						 temperature *tempp)
{
  //set the max altitude by finding the density at which the exosphere = n_species_min
  doubReal rmaxx = exosphere.r(n_species_min);
  
  setup_rmax_nCO2exo(rminn,
		     rexoo,
		     rmaxx,
		     rmindiffusionn,
		     n_species_exoo, // a good number is 10^5-6
		     nCO2exoo, //a good number is 10^9 (?)
		     tempp);
}

void thermosphere_exosphere::setup_rmax_nCO2rmin(doubReal rminn,
						 doubReal rexoo,
						 doubReal rmaxx,
						 doubReal rmindiffusionn,
						 doubReal n_species_exoo, 
						 doubReal nCO2rmin, //a good number is 2.6e13 (Chaufray 2008)
						 temperature *tempp)
{
  setup_rmax_nCO2exo(rminn,
		     rexoo,
		     rmaxx,
		     rmindiffusionn,
		     n_species_exoo, // a good number is 10^5-6
		     species_thermosphere->get_CO2_exobase_density(nCO2rmin, rminn, rexoo, tempp),
		     tempp);
}

void thermosphere_exosphere::setup_rmax_nCO2exo(doubReal rminn,
						doubReal rexoo,
						doubReal rmaxx,
						doubReal rmindiffusionn,
						doubReal n_species_exoo, 
						doubReal nCO2exoo, 
						temperature *tempp)
{
  rmin = rminn;
  rexo = rexoo;
  rmax = rmaxx;
    
  n_species_exo = n_species_exoo;
  nCO2exo = nCO2exoo;
  rmindiffusion = rmindiffusionn;

  temp = tempp;

  exosphere = chamberlain_exosphere(rexoo, temp->T_exo, n_species_exo, species_thermosphere->mass);
  species_thermosphere->get_thermosphere_density_arrays(n_species_exo,
							nCO2exo,
							rexo,
							rmin,
							
							exosphere.escape_flux,
							temp,
							
							log_nCO2_thermosphere,
							log_n_species_thermosphere,
							r_thermosphere,
							n_thermosphere_steps,
							/*get_interpolation_points = */true);

  CO2_exosphere = chamberlain_exosphere(rexoo, temp->T_exo, nCO2exo, mCO2);

#ifndef DNDEBUG
  //check the thermosphere values for any negatives
  for (unsigned int i=0; i < log_nCO2_thermosphere.size(); i++) {
    assert(exp(log_nCO2_thermosphere[i]) >= 0 && "densities must be positive.");
    assert(exp(log_n_species_thermosphere[i]) >= 0 && "densities must be positive.");
    assert(r_thermosphere[i] > 0 && "radii must be positive.");
  }
#endif
  
  //interpolate the densities in the thermosphere
  log_nCO2_thermosphere_spline = Linear_interp<doubReal>(r_thermosphere,log_nCO2_thermosphere);
  invlog_nCO2_thermosphere = Linear_interp<doubReal>(log_nCO2_thermosphere,r_thermosphere);

  log_n_species_thermosphere_spline = Linear_interp<doubReal>(r_thermosphere,log_n_species_thermosphere);
  invlog_n_species_thermosphere = Linear_interp<doubReal>(log_n_species_thermosphere,r_thermosphere);

  n_species_rmindiffusion = n_species(rmindiffusion);
  nCO2rmindiffusion = nCO2(rmindiffusion);
    
  //now get the interpolation points in the exosphere
  exosphere_step_logr = (log(rmax) - log(rexo))/(n_exosphere_steps - 1.);
  for (int iexo = 0; iexo < n_exosphere_steps; iexo++) {
    log_r_exosphere.push_back( log(rexo) + iexo * exosphere_step_logr );
    assert(log_r_exosphere.back() > 0 && "radii must be positive");
    log_n_species_exosphere.push_back( log( exosphere( exp( log_r_exosphere[iexo] ) ) ) );
    assert(exp(log_n_species_exosphere.back()) > 0 && "densities must be positive");
  }
  log_n_species_exosphere_spline = Linear_interp<doubReal>(log_r_exosphere,log_n_species_exosphere);
  invlog_n_species_exosphere = Linear_interp<doubReal>(log_n_species_exosphere,log_r_exosphere);


  exosphere_step_CO2_logr = (log(CO2_exo_zero_level) - log(rexo))/(n_exosphere_steps - 1.);
  for (int iexo = 0; iexo < n_exosphere_steps; iexo++) {
    log_r_CO2_exosphere.push_back( log(rexo) + iexo * exosphere_step_CO2_logr );
    assert(log_r_CO2_exosphere.back() > 0 && "radii must be positive");
    log_n_CO2_exosphere.push_back( log( CO2_exosphere( exp( log_r_CO2_exosphere[iexo] ) ) ) );
    assert(exp(log_n_CO2_exosphere.back()) > 0 && "densities must be positive");
  }
  log_n_CO2_exosphere_spline = Linear_interp<doubReal>(log_r_CO2_exosphere,log_n_CO2_exosphere);
  invlog_n_CO2_exosphere = Linear_interp<doubReal>(log_n_CO2_exosphere,log_r_CO2_exosphere);


  init=true;
}


doubReal thermosphere_exosphere::nCO2(const doubReal &r) const {
  if (r>CO2_exo_zero_level)
    return 0.0;
  else if (r>rexo && CO2_exo_zero_level!=rexo)
    return exp(log_n_CO2_exosphere_spline(log(r)));
  else {
    assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");
    return exp(log_nCO2_thermosphere_spline(r));
  }
}
doubReal thermosphere_exosphere::n_absorber(const doubReal &r) const  {
  return nCO2(r);
}

doubReal thermosphere_exosphere::n_species(const doubReal &r) const {
  if (r>=rexo)
    return exp(log_n_species_exosphere_spline(log(r)));
  else {
    if (r>=rmindiffusion)
      return exp(log_n_species_thermosphere_spline(r));
    else {
      assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");
      return n_species_rmindiffusion/nCO2rmindiffusion*nCO2(r);
    }
  }
}

doubReal thermosphere_exosphere::Temp(const doubReal &r) const {
  return temp->T(r);
}

doubReal thermosphere_exosphere::r_from_n_species(const doubReal &nsptarget) const {
  if (nsptarget==n_species_exo) {
    return rexo;
  } else if (nsptarget<n_species_exo) {
    return exp(invlog_n_species_exosphere(log(nsptarget)));
  } else if (nsptarget>n_species_rmindiffusion) {
    return invlog_nCO2_thermosphere(log(nsptarget*nCO2rmindiffusion/n_species_rmindiffusion));
  } else {
    return invlog_n_species_thermosphere(log(nsptarget));
  }
}

doubReal thermosphere_exosphere::nCO2_exact(const doubReal &r) const {
  if (r>rexo)
    return CO2_exosphere(r);
  else {
    assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");

    doubReal nCO2_tmp, n_species_tmp;
    species_thermosphere->get_thermosphere_density_exact(n_species_exo,
							 nCO2exo,
							 rexo,
							 
							 exosphere.escape_flux,
							 temp,
							 
							 r, 
							 nCO2_tmp,
							 n_species_tmp);
    return nCO2_tmp;
  }
}

doubReal thermosphere_exosphere::n_species_exact(const doubReal &r) const {
  if (r>=rexo)
    return exosphere(r);
  else {
    assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");
    if (r>=rmindiffusion) {

      doubReal nCO2_tmp, n_species_tmp;
      species_thermosphere->get_thermosphere_density_exact(n_species_exo,
							   nCO2exo,
							   rexo,
							   
							   exosphere.escape_flux,
							   temp,
							   
							   r, 
							   nCO2_tmp,
							   n_species_tmp);
      
      return n_species_tmp;
    } else 
      return n_species_exact(rmindiffusion)/nCO2_exact(rmindiffusion)*nCO2_exact(r);
  }
}



void thermosphere_exosphere::write_vector(std::ofstream &file, const std::string &preamble,
					  const vector<doubReal> &data) const {
  VectorXd write_out = Eigen::Map<const VectorXd>(data.data(),
						  data.size());

  file << preamble << write_out.transpose() << "\n";
}
  

void thermosphere_exosphere::save(std::string fname) const {

  std::ofstream file(fname.c_str());
  if (file.is_open())
    {

      file << "chaffin atmosphere for:\n"
	   << "         rexo = " << rexo << " cm,\n"
	   << "         Texo = " << temp->T_exo << " K,\n"
	   << "n_species_exo = " << n_species_exo << " cm-3,\n"
	   << "      nCO2exo = " << nCO2exo << " cm-3,\n\n";

      file << "Thermosphere is defined by solution to Krasnopolsky (2002) differential equation:\n";
      file << "Thermosphere interpolation is log-linear:\n";
      write_vector(file, "r [cm] = ", r_thermosphere);
      write_vector(file, "log(n_species) [cm-3] = ", log_n_species_thermosphere);
      write_vector(file, "log(nCO2) [cm-3] = ", log_nCO2_thermosphere);
      file << std::endl;

      file << "Exosphere is spherically symmetric Chamberlain (nCO2 assumed zero):\n";
      file << "Exosphere interpolation is log-log:\n";
      write_vector(file, "logr [cm] = ", log_r_exosphere);
      write_vector(file, "log(n_species) [cm-3] = ", log_n_species_exosphere);

      file.close();
    }
}
