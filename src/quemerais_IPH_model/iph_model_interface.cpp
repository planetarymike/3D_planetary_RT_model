#include "constants.hpp"
#include "iph_model_interface.hpp"
#include <cmath>

using std::vector;
using std::cos;
using std::sin;

//fortran function prototype
extern "C" {
  void background(const char* fname, const int* fname_length,
		  const float *lc,
		  const float *x_pos, const float *y_pos, const float *z_pos,
		  const int *n_los,
		  const float *x_look, const float *y_look, const float *z_look,
		  float *iph_b);
}

vector<Real> quemerais_iph_model(const string sfn_fname, // fully-qualified filename of Quemerais code input file
				 const Real &g_lya, //Lyman alpha g factor at Mars
				 const std::vector<Real> &marspos, //position of Mars in ecliptic coordinates [AU]
				 const vector<Real> &ra, const vector<Real> &dec) {
  //this routine computes the ecliptic look direction corresponding to
  //a given RA, dec, and feeds the appropriate parameters to the
  //Quemerais code, which returns a computed brightness.

  Real Fsun = g_lya / lyman_alpha_cross_section_total;//now photons/s/cm2/Hz
  Fsun *= (marspos[0]*marspos[0] + marspos[1]*marspos[1] + marspos[2]*marspos[2]);//now at position of Earth
  Fsun *= clight/lyman_alpha_lambda/lyman_alpha_lambda/1e8; //not photons/s/cm2/A
  
  int n_los = ra.size();
  vector<Real> iph_brightness;
  iph_brightness.resize(n_los);

  //fortran swap variables
  float lc_=Fsun;
  float x_pos_=marspos[0];
  float y_pos_=marspos[1];
  float z_pos_=marspos[2];
  float *x_look_ = new float[n_los];
  float *y_look_ = new float[n_los];
  float *z_look_ = new float[n_los];
  float *iphb_ = new float[n_los];

  for (int i_los=0; i_los<n_los; i_los++) {
    Real thisdec = M_PI/180*dec[i_los];
    Real thisra  = M_PI/180*ra[i_los];

    vector<Real> j2000vec(3);
    j2000vec[0]=cos(thisdec)*cos(thisra);
    j2000vec[1]=cos(thisdec)*sin(thisra);
    j2000vec[2]=sin(thisdec);

    //ecliptic coordinates are related to J2000 coordinates via a
    //negative rotation about the x-axis equal to the obliquity of the
    //Earth.
    Real eob=M_PI/180.*23.44;
    x_look_[i_los]=j2000vec[0];
    y_look_[i_los]=j2000vec[1]*cos(-eob)-j2000vec[2]*sin(-eob);
    z_look_[i_los]=j2000vec[2]*cos(-eob)+j2000vec[1]*sin(-eob);
  }

  // call the fortran code
  int sfn_fname_length_ = sfn_fname.length();
  background(sfn_fname.c_str(), &sfn_fname_length_,
	     &lc_,
	     &x_pos_,&y_pos_,&z_pos_,
	     &n_los,
	     &x_look_[0], &y_look_[0], &z_look_[0],
	     &iphb_[0]);

  for (int i_los=0; i_los<n_los; i_los++) {
    iph_brightness[i_los] = iphb_[i_los]/1000.; //iphb is in R, convert to kR
  }

  delete [] x_look_;
  delete [] y_look_;
  delete [] z_look_;
  delete [] iphb_;
  
  return iph_brightness;
}
