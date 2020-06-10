//csv_atmosphere.cpp -- read data from CSV files into splines for atmospheric parameters
#include "csv_atmosphere.hpp"

csv_atmosphere::csv_atmosphere(Real rminn, Real rexoo, Real rmaxx)
  : atmosphere(rminn, rexoo, rmaxx) { }

#include <fstream>
using std::string;

void csv_atmosphere::load_temperature(string fname) {
  Temp_alt.clear();
  Temp.clear();

  std::fstream file;
  file.open(fname.c_str());

  if(file.is_open()) {
    
    string line;
    string element;

    //discard header
    std::getline(file,line);
    
    while(std::getline(file,line)) {
      std::stringstream linestream(line);

      std::getline(linestream, element, ',');
      Temp_alt.push_back((Real) std::stod(element));

      std::getline(linestream, element, ',');
      Temp.push_back((Real) std::stod(element));
    }
  }

  T_spline = Linear_interp(Temp_alt,Temp);
}

void csv_atmosphere::load_densities(string fname) {
  n_alt.clear();
  lognCO2.clear();
  lognH.clear();

  std::fstream file;
  file.open(fname.c_str());

  if(file.is_open()) {
    
    string line;
    string element;

    //discard header
    std::getline(file,line);
    
    while(std::getline(file,line)) {
      std::stringstream linestream(line);

      std::getline(linestream, element, ',');
      n_alt.push_back((Real) std::stod(element));

      std::getline(linestream, element, ',');
      lognH.push_back((Real) log(std::pow(10.0,std::stod(element))));

      std::getline(linestream, element, ',');
      lognCO2.push_back((Real) log(std::pow(10.0,std::stod(element))));
    }
  }

  lognH_spline = Linear_interp(n_alt,lognH);
  inv_lognH_spline = Linear_interp(lognH,n_alt);
  lognCO2_spline = Linear_interp(n_alt,lognCO2);
}

void csv_atmosphere::scale_H(Real scale) {
  for (int i = 0; i < lognH.size(); i++)
    lognH[i]+=std::log(scale);
  lognH_spline = Linear_interp(n_alt,lognH);
  inv_lognH_spline = Linear_interp(lognH,n_alt);
}


Real csv_atmosphere::get_T(const Real &r) {
  return T_spline((r-rMars)/1e5);
}

Real csv_atmosphere::nCO2(const Real &r) {
  return exp(lognCO2_spline((r-rMars)/1e5));
}
Real csv_atmosphere::nCO2(const atmo_point pt) {
  return nCO2(pt.r);
}

Real csv_atmosphere::nH(const Real &r) {
  return exp(lognH_spline((r-rMars)/1e5));
}
Real csv_atmosphere::nH(const atmo_point pt) {
  return nH(pt.r);
}

Real csv_atmosphere::r_from_nH(Real nHtarget) {
  return inv_lognH_spline(log(nHtarget))*1e5 + rMars;
}

Real csv_atmosphere::sH_lya(const atmo_point pt) {
  assert(false && "not implemented in this class");
  return 0;
}

Real csv_atmosphere::sCO2_lya(const atmo_point pt) {
  assert(false && "not implemented in this class");
  return 0;
}
