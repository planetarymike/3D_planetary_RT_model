//atmo_vec.h -- basic definitions for points, rays, and vectors

#ifndef __atmo_vec_H
#define __atmo_vec_H

#include <vector>
#include <array>
#include <cmath>
#include <cassert>

using std::vector;
using std::array;

struct atmo_point {
  double x, y, z;
  double r, t, p;
  int i_voxel;//voxel index to which this point belongs
  bool init;

  atmo_point() {
    init=false;
    i_voxel = -1;
  }

  void rtp(const double rr, const double tt, const double pp) {
    r = rr; t = tt; p = pp;
    x = r*sin(t)*cos(p);
    y = r*sin(t)*sin(p);
    z = r*cos(t);
    init=true;
  }
  
  void xyz(const double xx, const double yy, const double zz) {
    x = xx; y = yy; z = zz;
    r=sqrt(x*x+y*y+z*z);
    t=acos(z/r);
    p=atan2(y,x);
    if (p<0)
      p+=2*pi;//puts the negative values on the right branch [0,2pi]
    init=true;
  }

  void set_voxel_index(const int ii) {
    assert(init && "point must be initialized to assign a voxel index.");
    i_voxel = ii;
  }  
};

struct atmo_ray {
  double t, p; // angles measured from the local vertical of the point in question
  double cost, sint; //needed for sphere intersections later on, might as well get them now
  int i_ray; //ray index, if the ray comes from the grid
  double domega; //solid angle belonging to this ray, if on grid; if not on a grid, 1.0 

  bool init;

  atmo_ray() {
    init=false;
    i_ray=-1;
    domega=1.0;
  }

  void tp(double tt, double pp) {
    t=tt;
    p=pp;

    cost = std::cos(t);
    sint = std::sin(t);
    init = true;
  }

  void set_ray_index(const int ii, double twt, double pwt) {
    if (init) {
      i_ray = ii;
      domega = sint*twt*pwt/4/pi;
    }
  }

  
};
  
struct atmo_vector {
  atmo_point pt;
  atmo_ray ray;
  
  double line_x, line_y, line_z;//cartesian vector elements 
  bool init;

  atmo_vector() {init=false;}

  atmo_vector(atmo_point ptt, atmo_ray rayy) : pt(ptt), ray(rayy)
  {
    //construct the cartesian elements from the point and ray
    //spherical coordinates.

    //Product of four rotation matrices.
    //probably worth refactoring to make this obvious
    line_x = cos(pt.p)*cos(ray.p)*cos(pt.t)*ray.sint
      +cos(pt.p)*ray.cost*sin(pt.t)
      -ray.sint*sin(ray.p)*sin(pt.p);
    line_y = ray.cost*sin(pt.p)*sin(pt.t)
      +ray.sint*cos(ray.p)*sin(pt.p)*cos(pt.t)
      +ray.sint*sin(ray.p)*cos(pt.p);
    line_z = ray.cost*cos(pt.t)-cos(ray.p)*ray.sint*sin(pt.t);

    init = true;
  }

  atmo_vector(atmo_point ptt, vector<double> vec) : atmo_vector(ptt,vec[0],vec[1],vec[2]) { }
  
  atmo_vector(atmo_point ptt, double line_xx, double line_yy, double line_zz) :
    pt(ptt) 
  {
    //this overload exists to construct rays from points toward the
    //sun, or rays from the spacecraft through the atmosphere. In both
    //cases the integration is along this line only so we can set
    //domega = 1. We need only set the theta values since these are
    //used in computing intersections.
    
    double line_mag=sqrt(line_xx*line_xx + line_yy*line_yy + line_zz*line_zz);
    line_x = line_xx/line_mag;
    line_y = line_yy/line_mag;
    line_z = line_zz/line_mag;

    ray.cost=(line_x*pt.x + line_y*pt.y + line_z*pt.z)/(pt.r);
    ray.sint=sqrt(1.0-ray.cost*ray.cost);

    ray.t=acos(ray.cost);
    ray.p=-1.0;
    
    ray.domega=1.0;
    ray.init=true;

    init=true;
  }

  atmo_point extend(double dist) {
    atmo_point retpt;

    retpt.xyz(pt.x+line_x*dist,
	      pt.y+line_y*dist,
	      pt.z+line_z*dist);

    return retpt;

      
  }
  
};

#endif

