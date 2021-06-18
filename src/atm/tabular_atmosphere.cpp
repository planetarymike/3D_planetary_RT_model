//tabular_atmosphere.cpp -- conpute densities from a table of atmospheric parameters
#include "tabular_atmosphere.hpp"

tabular_atmosphere::tabular_atmosphere(const double rminn, const double rexoo, const double rmaxx, 
				       const bool compute_exospheree/* = false*/,
				       const double m_speciess /* = mH */)
  : atmosphere(rminn,rexoo,rmaxx), compute_exosphere(compute_exospheree), m_species(m_speciess)
{ }

void tabular_atmosphere::load_log_species_density(const vector<double> &alt,
						  const vector<double> &log_n_species) {
  //log_n_species better be monotonic!
  bool monotonic = true;
  bool ascnd = log_n_species[1] > log_n_species[0] ? true : false;

  for (unsigned int i=1;i<log_n_species.size();i++)
    monotonic = monotonic && (ascnd ?
			      log_n_species[i] > log_n_species[i-1] :
			      log_n_species[i] < log_n_species[i-1]);
  assert(monotonic && "log_n_species must be monotonic and invertible");
  
  log_n_species_spline = Linear_interp<double>(alt,log_n_species);
  inv_log_n_species_spline = Linear_interp<double>(log_n_species,alt);

  check_init();
}

void tabular_atmosphere::load_log_absorber_density(const vector<double> &alt,
						   const vector<double> &log_n_absorber) {
  log_n_absorber_spline = Linear_interp<double>(alt,log_n_absorber);
  check_init();
}

void tabular_atmosphere::load_temperature(const vector<double> &alt,
					  const vector<double> &temp) {
  Temp_spline = Linear_interp<double>(alt,temp);

  check_init();
}

void tabular_atmosphere::check_init() {
  if (compute_exosphere && log_n_species_spline.n > 0 && Temp_spline.n > 0)
    init_exosphere();
  if (log_n_species_spline.n > 0 && Temp_spline.n > 0 && log_n_absorber_spline.n > 0)
    init = true;
}

void tabular_atmosphere::init_exosphere() {
    double n_species_exo = n_species(rexo);
    double T_exo = Temp(rexo);
    
    //set up the exosphere from the values at the exobase
    exosphere = chamberlain_exosphere(rexo, T_exo, n_species_exo, m_species);
};

double tabular_atmosphere::n_species(const double &r) const {
  assert(log_n_species_spline.n > 0 && "n_species must be initialized!");

  if (compute_exosphere && r>rexo) {
    return exosphere.n(r);
  } else  
    return exp(log_n_species_spline((r-rMars)/1e5));
}

double tabular_atmosphere::r_from_n_species(const double &n_species_target) const {
  assert(inv_log_n_species_spline.n > 0 && "n_species must be initialized!");
  if (compute_exosphere && n_species_target < n_species(rexo))
    return exosphere.r(n_species_target);
  else
    return inv_log_n_species_spline(log(n_species_target))*1e5 + rMars;
}

double tabular_atmosphere::Temp(const double &r) const {
  assert(Temp_spline.n > 0 && "Temp must be initialized!");
  if (compute_exosphere && r>rexo)
    return Temp(rexo);
  else
    return Temp_spline((r-rMars)/1e5);
}

double tabular_atmosphere::n_absorber(const double &r) const {
  assert(log_n_absorber_spline.n > 0 && "n_absorber must be initialized!");
  if (compute_exosphere && r>rexo)
    return 0.0;
  else
    return exp(log_n_absorber_spline((r-rMars)/1e5));
}
