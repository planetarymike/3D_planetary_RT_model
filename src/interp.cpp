#include "interp.hpp"

int Base_interp::locate(const Real x)
// given a value x, return a value j such that x is (insofar as
// possible) centered in the subrange xx[j..j+mmm-1], where xx is the
// stored pointer. The values in xx must be monotonic, either
// increasing or decreasing. The returned value is not less than 0,
// nor greater than n-1).
{
  int ju, jm, jl;
  //  std::cout << "n = " << n << ", mmm = " << mmm << "\n";
  if (n < 2 || mmm < 2 || mmm > n) {
    //    std::cout << "n = " << n << ", mmm = " << mmm << std::endl;
    assert(false && "locate size error");
  }

  bool ascnd = (xx[n-1] >= xx[0]); // true if ascending order of
				   // table, false otherwise
  jl = 0;//lower and upper size limits
  ju = n-1;
  while (ju-jl > 1) {  // until the appropriate range is found
    jm = (ju+jl) >> 1; // find the midpoint
    if ((x >= xx[jm]) == ascnd)
      jl = jm; // and replace either the lower limit
    else
      ju = jm; // or the upper limit, as appropriate
  }
  cor = abs(jl-jsav) > dj ? 0 : 1; // decide whether to use hunt() or
				   // locate() next time
  jsav = jl;
  return std::max(0,std::min(n-mmm,jl-((mmm-2)>>1)));
}


int Base_interp::hunt(const Real x)
// given a value x, return a value j such that x is (insofar as
// possible) centered in the subrange xx[j..j+mmm-1], where xx is the
// stored pointer. The values in xx must be monotonic, either
// increasing or decreasing. The returned value is not less than 0,
// nor greater than n-1).
{
  int jl = jsav, jm, ju, inc = 1;
  //  std::cout << "n = " << n << ", mmm = " << mmm << std::endl;
  if (n < 2 || mmm < 2 || mmm > n) {
    //    std::cout << "n = " << n << ", mmm = " << mmm << std::endl;
    assert(false && "hunt size error");
  }
  bool ascnd = (xx[n-1] >= xx[0]); // does table ascend?
  if (jl < 0 || jl > n-1) { // input guess not useful. go directly to
			    // bisection
    jl = 0;
    ju = n-1;
  } else {
    if ((x >= xx[jl]) == ascnd) { // hunt up:
      for(;;) {
	ju = jl + inc;
	if (ju >= n-1) { ju = n-1; break; } // off end of table
	else if ((x < xx[ju]) == ascnd) break; // found bracket
	else {
	  jl = ju;
	  inc += inc;
	}
      }
    } else { // hunt down:
      ju = jl;
      for (;;) {
	jl = jl - inc;
	if (jl <= 0) {jl = 0; break;} // off end of table
	else if ((x >= xx[jl]) == ascnd) break; // found bracket
	else {
	  ju = jl;
	  inc += inc;
	}
      }
    }
  }
  while (ju - jl > 1) { // hunt is done, so begin final bisection
    jm = (ju+jl) >> 1;
    if ((x >= xx[jm]) == ascnd)
      jl = jm;
    else
      ju = jm;
  }
  cor = abs(jl - jsav) > dj ? 0 : 1; // decide whether to use hunt or
				     // locate next time
  jsav = jl;
  return std::max(0,std::min(n-mmm,jl-((mmm-2)>>1)));
}
