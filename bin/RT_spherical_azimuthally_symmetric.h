// RT_spherical_azimuthally_symmetric.h -- spherical grid with symmetry about the Mars-Sun line

#ifndef __RT_spherical_azimuthally_symmetric
#define __RT_spherical_azimuthally_symmetric

#include <RT_grid.h>
#include "coordinate_generation.h"
#include "boundaries.h"
#include "atmosphere.h"
#include "interp.h"
#include "intersections.h"
#include <type_traits> 
#include <fstream>
#include <cmath>
#include <string>

struct spherical_azimuthally_symmetric_grid : RT_grid
{

  const int r_dimension;
  int n_radial_boundaries;
  int rmethod;
  static const int rmethod_altitude = 0;
  static const int rmethod_lognH = 1;
  
  vector<double> radial_boundaries;
  vector<double> pts_radii;
  vector<sphere> radial_boundary_spheres;

  const int sza_dimension;
  int n_sza_boundaries;
  int szamethod;
  static const int szamethod_uniform = 0;
  static const int szamethod_uniform_cos = 1;

  vector<double> sza_boundaries;
  vector<double> pts_sza;
  vector<cone> sza_boundary_cones;


  vector<int> n_boundaries;
  intersection_writer saver;
  bool save_intersections;

  int n_theta;
  vector<double> ray_theta;
  int n_phi;
  vector<double> ray_phi;

  
 spherical_azimuthally_symmetric_grid(int n_emissions,
				      influence &transmissionn)
   : RT_grid(n_emissions,transmissionn),
    r_dimension(0), sza_dimension(1)
    {
      n_dimensions = 2;
      n_boundaries.resize(n_dimensions);
      sun_direction = {0.,0.,1.};
      rmethod = rmethod_altitude;
      szamethod = szamethod_uniform;
      
      save_intersections = false;
      saver.fname = "intersections.dat";
    }

  void setup_voxels(int n_radial_boundariess,
		    int n_sza_boundariess,
		    atmosphere &atm)    
  {
    n_radial_boundaries = n_radial_boundariess;
    
    rmin = atm.rmin;
    rmax = atm.rmax;

    assert((rmethod == rmethod_altitude || rmethod == rmethod_lognH)
	   && "rmethod must match a defined radial points method");

    // don't define a tau radial points method; tau < 0.1 is
    // important and max(tau) > 10; this leads to many required
    // gridpoints
    
    if (rmethod == rmethod_altitude)
      get_radial_log_linear_points(radial_boundaries, n_radial_boundaries,
				   atm.rmin, atm.rexo, atm.rmax);
    if (rmethod == rmethod_lognH) {
      double lognH_max = log(atm.nH(atm.rmin));
      double lognH_min = log(atm.nH(atm.rmax));
      double lognH_step = (lognH_max-lognH_min)/(n_radial_boundaries-1.);
      
      for(int i=0;i<n_radial_boundaries;i++) {
	double nH_target=exp(lognH_max-i*lognH_step);
	radial_boundaries.push_back(atm.r_from_nH(nH_target));
      }
    }
    
    pts_radii.resize(n_radial_boundaries-1);
    for (int i=0; i<n_radial_boundaries-1; i++) {
      pts_radii[i]=sqrt(radial_boundaries[i]*radial_boundaries[i+1]);
    }

    radial_boundary_spheres.resize(n_radial_boundaries);
    for (int i=0; i<n_radial_boundaries; i++) 
      radial_boundary_spheres[i].set_radius(radial_boundaries[i]);



    n_sza_boundaries = n_sza_boundariess;

    assert((szamethod == szamethod_uniform || szamethod == szamethod_uniform_cos)
	   && "szamethod must match a defined sza points method");
    
    
    if (szamethod == szamethod_uniform) {
      double sza_spacing = pi / (n_sza_boundaries - 2.);
      sza_boundaries.resize(n_sza_boundaries);
      for (int i=0;i<n_sza_boundaries;i++) {
	sza_boundaries[i]=(i-0.5)*sza_spacing;
      }
    }
    if (szamethod == szamethod_uniform_cos) {
      double cos_sza_spacing = 2.0 / (n_sza_boundaries - 2.);
      using std::acos;
      
      sza_boundaries.resize(n_sza_boundaries);
      sza_boundaries[0] = -acos(1.0-0.5*cos_sza_spacing);
      for (int i=1;i<n_sza_boundaries-1;i++) {
	sza_boundaries[i]=acos(1.0-(i-0.5)*cos_sza_spacing);
      }
      sza_boundaries[n_sza_boundaries-1] = pi + acos(1.0-0.5*cos_sza_spacing);
    }
    

    pts_sza.resize(n_sza_boundaries-1);
    for (unsigned int i=0;i<pts_sza.size();i++) {
      pts_sza[i]=0.5*(sza_boundaries[i] + sza_boundaries[i+1]);
    }
    
    sza_boundary_cones.resize(n_sza_boundaries-2);
    for (int i=0;i<n_sza_boundaries-2;i++)
      sza_boundary_cones[i].set_angle(sza_boundaries[i+1]);



    int n_voxelss = pts_radii.size()*pts_sza.size();
    n_boundaries[r_dimension] = n_radial_boundaries-1;
    n_boundaries[sza_dimension] = n_sza_boundaries-1;
    init_arrays(n_voxelss);
    
    pts.resize(n_voxels);
    int ivoxel;
    for (unsigned int i=0; i<pts_radii.size(); i++) {
      for (unsigned int j=0;j<pts_sza.size();j++) {
	ivoxel = i*(n_sza_boundaries-1) + j;
	
	pts[ivoxel].rtp(pts_radii[i],pts_sza[j],0.);
	pts[ivoxel].set_voxel_index(ivoxel);
      }
    }

    if (save_intersections)
      saver.save_coordinates(radial_boundaries,sza_boundaries);
    
    voxels_init=true;
    if (rays_init)
      grid_init=true;
  }
  
  void setup_rays(int n_thetaa, int n_phii)
  {
    n_theta = n_thetaa;
    n_phi = n_phii;

    vector<double> ray_weights_theta;
    gauss_quadrature_points(ray_theta,ray_weights_theta,0,pi,n_theta);
    
    double phi_spacing = 2*pi/n_phi;
    ray_phi.resize(n_phi);
    for (int i=0;i<n_phi;i++)
      ray_phi[i] = (i+0.5)*phi_spacing;
    
    n_rays = n_theta*n_phi;
    rays.resize(n_rays);
    int iray;
    for (int i=0;i<n_theta;i++) {
      for (int j=0;j<n_phi;j++) {
	iray = i * n_phi + j;
	rays[iray].tp(ray_theta[i],ray_phi[j]);
	rays[iray].set_ray_index(iray, ray_weights_theta[i], phi_spacing);
      }
    }
    
    rays_init=true;
    if (voxels_init)
      grid_init = true;
  }
  
  int indices_to_voxel(const vector<int> &indices) {
    unsigned int r_idx = indices[r_dimension];
    unsigned int sza_idx = indices[sza_dimension];
    
    if ((r_idx   < 0) || (  r_idx > pts_radii.size()-1) ||
	(sza_idx < 0) || (sza_idx > pts_sza.size()-1))
      return -1;
    else
      return indices[r_dimension]*(n_sza_boundaries-1)+indices[sza_dimension];
  }

  vector<int> voxel_to_indices(const int i_voxel) {
    vector<int> v(n_dimensions,-1);
    if ((i_voxel < 0) || (i_voxel > n_voxels-1)) {
      v[r_dimension]=-1;
      v[sza_dimension]=-1;
    } else {
      v[r_dimension]=i_voxel / (n_sza_boundaries-1);
      v[sza_dimension]=i_voxel % (n_sza_boundaries-1);
    }
    
    return v;
  }

  int find_coordinate_index(double pt_coord, vector<double> boundaries) {
    unsigned int i;
    
    for (i=0;i<boundaries.size();i++)
      if (pt_coord < boundaries[i])
	break;

    i--;

    assert((boundaries.back()<=pt_coord ||
	    pt_coord < boundaries[0] ||
	    (boundaries[i]<=pt_coord &&pt_coord<boundaries[i+1]))
	   && "we have found the appropriate point index");
    
    return i;
  }

  vector<int> point_to_indices(const atmo_point pt) {
    vector<int> indices(n_dimensions,-1);
    indices[r_dimension] = find_coordinate_index(pt.r, radial_boundaries);
    indices[sza_dimension] = find_coordinate_index(pt.t, sza_boundaries);
    return indices;
  }
  

  boundary_intersection_stepper ray_voxel_intersections(const atmo_vector &vec) {
    
    boundary_set boundaries(n_dimensions);

    //define the origin
    boundary origin(n_dimensions);
    if (vec.pt.i_voxel == -1) {
      origin.entering_indices = point_to_indices(vec.pt);
      origin.entering = indices_to_voxel(origin.entering_indices);
    } else {
      origin.entering = vec.pt.i_voxel;
      origin.entering_indices = voxel_to_indices(origin.entering);
    }
    origin.distance = 0.0;
    boundaries.append(origin);

    //do the intersections for each coordinate
    vector<double> temp_distances;
    for (unsigned int ir=0;ir<radial_boundary_spheres.size();ir++) {
      temp_distances = radial_boundary_spheres[ir].intersections(vec);
      if (temp_distances.size()>0) 
	boundaries.add_intersections(vec.pt.r, r_dimension,
				     ir, radial_boundaries[ir], temp_distances);
    }

    for (unsigned int isza=0;isza<sza_boundary_cones.size();isza++) {
      temp_distances = sza_boundary_cones[isza].intersections(vec);
      if (temp_distances.size() > 0)
      	boundaries.add_intersections(vec.pt.t, sza_dimension,
      				     isza+1, sza_boundaries[isza+1], temp_distances);
    }

    //sort the list of intersections by distance & trim
    boundaries.sort();
    boundaries.propagate_indices();
    boundaries.assign_voxel_indices(this, &spherical_azimuthally_symmetric_grid::indices_to_voxel);
    boundaries.trim();
    boundaries.check(n_boundaries,n_voxels);
    
    if (save_intersections)
      saver.append_intersections(vec,boundaries);

    return boundary_intersection_stepper(vec, boundaries);
  }

  VectorXd sza_slice(VectorXd quantity, int i_sza) {
    VectorXd ret;
    vector<int> indices(n_dimensions,i_sza);
    
    ret.resize(n_radial_boundaries-1);
    for (int i=0;i<n_radial_boundaries-1;i++) {
      indices[0]=i;
      ret(i) = quantity(indices_to_voxel(indices));
    }
    
    return ret;
  }


  
  class spherical_azimuthally_symmetric_interp : public qinterp {
    Bilinear_interp terp;
  public:
    spherical_azimuthally_symmetric_interp() : qinterp() { }

    template <class T>
    spherical_azimuthally_symmetric_interp(VectorXd &input,
					   T *parent)
    {
      MatrixXd data;
      vector<double> log_radii;
      vector<double> szas;

      int nrpts = parent->pts_radii.size();
      int nlogradii = nrpts+2;
      int nsza = parent->pts_sza.size();
      
      data.resize(nlogradii,nsza);

      log_radii.push_back(log(parent->radial_boundaries[0]));
      for (int i=0;i<nrpts;i++) 
	log_radii.push_back(log(parent->pts_radii[i]));
      log_radii.push_back(log(parent->radial_boundaries.back()));

      for (int j=0;j<nsza;j++) 
	szas.push_back(parent->pts_sza[j]);

      vector<int> idx(2,0);
      for (int i=0;i<nlogradii;i++) {
	for (int j=0;j<nsza;j++) {
	  if (i==0) {
	    idx[0]=0; idx[1]=j;
	  } else if (i==nlogradii-1) {
	    idx[0]=nrpts-1; idx[1]=j;
	  } else {
	    idx[0]=i-1; idx[1]=j;
	  }
	  data(i,j) = input(parent->indices_to_voxel(idx));
	}
      }

      terp = Bilinear_interp(log_radii, szas, data);
    }

    double interp(const atmo_point pt) { 
      return terp.interp(log(pt.r),pt.t);
    }

  };

  class spherical_azimuthally_symmetric_grid_interpolator {
  public:
    spherical_azimuthally_symmetric_interp dtau_species;
    spherical_azimuthally_symmetric_interp dtau_absorber;
    spherical_azimuthally_symmetric_interp abs;
    spherical_azimuthally_symmetric_interp sourcefn;
  }; 

  vector<spherical_azimuthally_symmetric_grid_interpolator> internal_interp;
  
  void setup_interp() {
    interp.resize(n_emissions);
    internal_interp.resize(n_emissions);

    for (int i_emission=0;i_emission<n_emissions;i_emission++) {
      internal_interp[i_emission].dtau_species  = spherical_azimuthally_symmetric_interp(dtau_species[i_emission],
											 this);
      interp[i_emission].dtau_species = &internal_interp[i_emission].dtau_species;
      
      internal_interp[i_emission].dtau_absorber = spherical_azimuthally_symmetric_interp(dtau_absorber[i_emission],
											 this);
      interp[i_emission].dtau_absorber = &internal_interp[i_emission].dtau_absorber;
      
      internal_interp[i_emission].abs           = spherical_azimuthally_symmetric_interp(abs[i_emission],
											 this);
      interp[i_emission].abs = &internal_interp[i_emission].abs;
      
      internal_interp[i_emission].sourcefn      = spherical_azimuthally_symmetric_interp(sourcefn[i_emission],
											 this);
      interp[i_emission].sourcefn = &internal_interp[i_emission].sourcefn;
    }
    interp_init=true;
  }
    
  void save_S(string fname) {
    std::ofstream file(fname.c_str());
    if (file.is_open())
      {

	VectorXd r_boundaries_write_out = Eigen::Map<VectorXd>(radial_boundaries.data(),
							       radial_boundaries.size());

	file << "radial boundaries [cm]: " << r_boundaries_write_out.transpose() << "\n\n";

	VectorXd r_pts_write_out = Eigen::Map<VectorXd>(pts_radii.data(),
							pts_radii.size());

	file << "pts radii [cm]: " << r_pts_write_out.transpose() << "\n\n";
       	
	VectorXd sza_boundaries_write_out = Eigen::Map<VectorXd>(sza_boundaries.data(),
								 sza_boundaries.size());

	file << "sza boundaries [rad]: " << sza_boundaries_write_out.transpose() << "\n\n";
	
	VectorXd sza_pts_write_out = Eigen::Map<VectorXd>(pts_sza.data(),
							  pts_sza.size());

	file << "pts sza [rad]: " << sza_pts_write_out.transpose() << "\n\n";
       	
	for (int i=0; i<n_emissions; i++) {
	  file << "For " << emission_names[i] << "\n";
	  for (unsigned int j=0; j<pts_sza.size(); j++) {
	    file << "  For SZA = " << pts_sza[j] << ": \n" 
		 << "    Species density [cm-3]: "
		 <<      sza_slice(species_density[i],j).transpose() << "\n"

		 << "    Species single scattering tau: " 
		 <<	 sza_slice(tau_species_single_scattering[i],j).transpose() << "\n"

		 << "    Species cross section [cm2]: " 
		 <<      sza_slice(species_sigma[i],j).transpose() << "\n"

		 << "    Absorber density [cm-3]: " 
		 <<      sza_slice(absorber_density[i],j).transpose() << "\n"

		 << "    Absorber single scattering tau: " 
		 <<	 sza_slice(tau_absorber_single_scattering[i],j).transpose() << "\n"

		 << "    Source function: " 
		 <<      sza_slice(sourcefn[i],j).transpose() << "\n\n";
	  }
	}
      }
    file.close();
  }

  
};




#endif
