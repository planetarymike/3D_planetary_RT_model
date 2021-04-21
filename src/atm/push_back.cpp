#include "push_back.hpp"

void push_back_quantities::operator()( const vector<double> &x , double r )
{
  for(unsigned int i=0; i<x.size(); i++)
    (*(vec_ptrs[i])).push_back(x[i]);
  (*(vec_ptrs.back())).push_back( r );
}
