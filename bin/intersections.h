//intersections.h -- routines for computing intersections between lines and geometric primitives

#ifndef __INTERSECTIONS_H
#define __INTERSECTIONS_H

#include "atmo_vec.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

using std::vector;

struct plane {
  double z;

  vector<double> intersections(const atmo_vector & vec) {
    vector<double> distances;

    if (vec.line_z != 0) {
      double d = (z-vec.pt.z)/vec.line_z;
      if (d > 0)
	distances.push_back(d);
    }

    return distances;    
  }
  
};



class sphere {
private:
  double r;
  double r2;

public:
  void set_radius(const double &rr) {
    r=rr;
    r2=r*r;
  }

  vector<double> intersections(const atmo_vector & vec) {
    vector<double> distances;

    // we solve the quadratic t^2 - 2*(vec.r * vec.cost)*t + (vec.r^2 - r^2) == 0
    double B = vec.pt.r * vec.ray.cost;
    double C = vec.pt.r*vec.pt.r - r2;

    double discr = B*B-C;

    //we can ignore cases where the ray misses or is tangent
    if (discr > 0) {
      discr = std::sqrt(discr);
      
      double d0 = -B - discr;
      if (d0>0) 
	distances.push_back(d0);

      double d1 = -B + discr;
      if (d1>0) 
	distances.push_back(d1);
    }
    return distances;
  }

  
};


class cone {
private:  
  double angle;
  double cosangle;
  double cosangle2;

  bool samesign(double a, double b) {
    if ((a>=0&&b>=0) || (a<0&&b<0))
      return true;
    else
      return false;
  }

  bool is_zero(double a) {
    double tol = 1e-10;
    if (a>tol || -a < -tol)
      return false;
    else
      return true;
  }
  
public:
  
  void set_angle(const double &a) {
    angle=a;
    cosangle=std::cos(angle);
    cosangle2=cosangle*cosangle;
  }

  vector<double> intersections(const atmo_vector & vec) {
    vector<double> distances;

    double A = vec.line_z * vec.line_z - cosangle2;
    double B = vec.pt.z * vec.line_z - vec.pt.r * vec.ray.cost * cosangle2;
    double C = vec.pt.z * vec.pt.z - vec.pt.r * vec.pt.r * cosangle2;
    
    if (!is_zero(A)) {
      double discr = B*B-C;

      if (discr > 0) {
	discr = std::sqrt(discr);
	
	double d0 = (-B - discr)/A;
	if (d0 > 0 && samesign(vec.pt.z + d0*vec.line_z, cosangle) ) 
	  distances.push_back(d0);

	double d1 = (-B + discr)/A;
	if (d1 > 0 && samesign(vec.pt.z + d1*vec.line_z, cosangle)) 
	  distances.push_back(d1);
      }
    } else {
      double d = -C/B;
      if (d>0 && samesign(vec.pt.z + d*vec.line_z, cosangle)) 
	distances.push_back(d);
    }
    
    std::sort(distances.begin(),distances.end());
    return distances;
  }
};


struct intersection_writer {
  string fname;
  std::ofstream file;

  intersection_writer() { }
  
  intersection_writer(string fnamee) : fname(fnamee) { }

  ~intersection_writer() {
    if (file.is_open()) {
      file.close();
    }
  }
  
  void save_coordinates(vector<double> radial_boundaries,
			vector<double> sza_boundaries) {
    file.open(fname.c_str());
    if (file.is_open()) {
      VectorXd r_boundaries_write_out = Eigen::Map<VectorXd>(radial_boundaries.data(),
							     radial_boundaries.size());
    
      file << "radial boundaries [cm]:\n" << r_boundaries_write_out.transpose() << "\n";
    
      VectorXd sza_boundaries_write_out = Eigen::Map<VectorXd>(sza_boundaries.data(),
							       sza_boundaries.size());
    
      file << "sza boundaries [rad]:\n" << sza_boundaries_write_out.transpose() << "\n\n";
    }
  }

  void append_intersections(atmo_vector vec,
			    boundary_set boundaries) {
    if (file.is_open()) {
      
      file << "ray origin (xyz) [cm]:\n"
	   << vec.pt.x << "  "
	   << vec.pt.y << "  "
	   << vec.pt.z << "\n";
      
      file << "ray direction (xyz):\n"
	   << vec.line_x << "  "
	   << vec.line_y << "  "
	   << vec.line_z << "\n";
      
      
      file << "Intersections [t, r_bound, sza_bound]\n";
      for (unsigned int i=0;i<boundaries.size();i++) {
	file << boundaries[i].distance << "  "
	     << boundaries[i].entering_indices[0] << "  "
	     << boundaries[i].entering_indices[1] << "\n";
      }
      file << "\n";
    }
  }
};
  


#endif
