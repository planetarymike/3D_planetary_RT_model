// Real.cpp

#include <cstdlib>
#include <iostream>
#include "Real.hpp"

CUDA_CALLABLE_MEMBER 
void assert_finite(const Real &value) {
#if !defined(NDEBUG)
  bool isnan = true;
  isnan *= std::isnan(value);
  
  if (isnan) {
    printf("value must be finite.\n");
    assert(false);
  }
#endif
}

CUDA_CALLABLE_MEMBER 
void assert_positive(const Real &value) {
#if !defined(NDEBUG)
  Real test_val = 1.0;
  test_val *= value;
  
  if (!(test_val >= REAL(0.))) {
    printf("value must be positive.\n");
    assert(false);
  }
#endif
}

CUDA_CALLABLE_MEMBER 
void assert_probability(const Real &value) {
#if !defined(NDEBUG)
  Real test_val = 1.0;
  test_val *= value;
  
  if (!(REAL(0.) <= test_val && test_val <= REAL(1.))) {
    printf("value must represent a probability.\n");
    assert(false);
  }
#endif
}

CUDA_CALLABLE_MEMBER 
void assert_leq(const Real &value, const Real &limit) {
#if !defined(NDEBUG)
  Real test_val = 1.0;
  test_val *= value;

  if (!(test_val >= 0 &&
	(test_val <= limit ||
	 std::abs(1.0-limit/test_val) < EPS))) {
    printf("value must be less than limit.\n");
    assert(false);
  }
#endif
}

CUDA_CALLABLE_MEMBER 
void assert_small_contribution(const Real &expected_small,
			       const Real comparison_value/*=1.0*/,
			       const Real abstol/*=1e-6*/,
			       const Real reltol/*=1e-2*/) {
#if !defined(NDEBUG)
  Real comp_val = 1.0;
  comp_val *= comparison_value;

  if (comp_val > abstol && !(reltol*comp_val > expected_small)) {
    printf("value must be small compared to reference.\n");
    assert(false);
  }
#endif
}

