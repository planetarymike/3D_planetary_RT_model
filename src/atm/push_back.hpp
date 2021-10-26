#ifndef __PUSH_BACK_H
#define __PUSH_BACK_H

#include "Real.hpp"
#include <vector>
using std::vector;

struct push_back_quantities
{
  vector< vector<doubReal>* > vec_ptrs;

  //this is copypasta that makes a variadic constructor that stores an
  //arbitrary number of passed vector pointers in vec_ptrs
  template<bool...> struct bool_pack{};
  template<class... Ts>
  using conjunction = std::is_same<bool_pack<true,Ts::value...>, bool_pack<Ts::value..., true>>;
  template<typename... Ts>
  using AllVecs = typename std::enable_if<conjunction<std::is_convertible<Ts, vector<doubReal>*>...>::value>::type;
  template<typename... Ts, typename = AllVecs<Ts...>>
  push_back_quantities(Ts... v0) {
    vector<vector<doubReal>*> vecs = { v0... };
    for (auto&& v : vecs) {
      vec_ptrs.push_back(v);
    }
  }

  void operator()( const vector<doubReal> &x , doubReal r );
};



#endif
