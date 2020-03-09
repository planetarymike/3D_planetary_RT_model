//geometric_tools_interface.h

#ifndef __GEOMETRIC_TOOLS_INTERFACE_H
#define __GEOMETRIC_TOOLS_INTERFACE_H

#include <Mathematics/IntrRay3Plane3.h>
#include <Mathematics/IntrRay3Sphere3.h>
#include <Mathematics/IntrRay3Cone3.h>
#include <cmath>

using gte::FIQuery;
using gte::Plane3;
using gte::Sphere3;
using gte::Cone3;
using gte::QFNumber;
using gte::Ray3;

double tol = 1e-5;


typedef FIQuery<double, Ray3<double>, Plane3<double>> gte_ray_plane_query;
class ray_plane_query : gte_ray_plane_query {
public:
  bool intersect;
  vector<double> distance;

  void operator()(Ray3<double> const& ray, Plane3<double> const& plane) {
    gte_ray_plane_query::Result Result = gte_ray_plane_query::operator()(ray,plane);
    intersect=false;
    distance.clear();

    if (Result.intersect) {
      assert(Result.numIntersections == 1
	     && "Only one intersection is possible between a ray ans a plane.");
      intersect = true;
      distance.resize(1);
      distance[0]=Result.parameter;
    }
  }
};

typedef FIQuery<double, Ray3<double>, Sphere3<double>> gte_ray_sphere_query;
class ray_sphere_query : gte_ray_sphere_query {
public:
  bool intersect;
  vector<double> distances;

  void operator()(Ray3<double> const& ray, Sphere3<double> const& sphere) {
    gte_ray_sphere_query::Result Result = gte_ray_sphere_query::operator()(ray,sphere);
    intersect=false;
    distances.clear();
    
    if (Result.intersect) {
      assert((Result.numIntersections==1||Result.numIntersection==2)
	     && "only 1 or 2 intersections are possible between a ray and a sphere");

      
      vector<double> test_distances;
      if (Result.numIntersections == 1)
	test_distances.push_back(Result.parameter[0]);
      if (Result.numIntersections == 2) {
	test_distances.push_back(Result.parameter[0]);
	test_distances.push_back(Result.parameter[1]);
      }
      
      for (auto d: test_distances) {
	if ((d>tol) && !(distances.size() > 0 && d-distances.back()>tol) )
	  distances.push_back(d);
      }
      if (distances.size() > 0)
	intersect=true;
    }
  }
};

typedef FIQuery<double, Ray3<double>, Cone3<double>> gte_ray_cone_query;
class ray_cone_query : gte_ray_cone_query {
public:
  bool intersect;
  vector<double> distances;

  double QF_to_double(QFNumber<double, 1> t) {
    if (t.x[1] != 0.)
      return t.x[0]+t.x[1]*std::sqrt(t.d);
    else
      return t.x[0];
  }
  
  void operator()(Ray3<double> const& ray, Cone3<double> const& cone) {
    gte_ray_cone_query::Result Result = gte_ray_cone_query::operator()(ray,cone);
    intersect=false;
    distances.clear();

    if (Result.intersect) {
      vector<double> test_distances;
	if (Result.type == Result.isSegment) {
	  test_distances.push_back(QF_to_double(Result.t[0]));
	  test_distances.push_back(QF_to_double(Result.t[1]));
	} else if (Result.type == Result.isRayPositive) {
	  test_distances.push_back(QF_to_double(Result.t[0]));
	} else if (Result.type == Result.isRayPositive) {
	  test_distances.push_back(QF_to_double(Result.t[1]));
	}
	
	for (auto d: test_distances) {
	  if ((d>tol) && !(distances.size() > 0 && d-distances.back()>tol) )
	    distances.push_back(d);
	}
	if (distances.size() > 0)
	  intersect=true;
    }
  }
};









#endif
