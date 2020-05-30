//my_clock.hpp

#ifndef __my_clock_H
#define __my_clock_H

#include "Real.hpp"
#include <ctime>
#include <iostream>
#include <string>

struct my_clock {
  clock_t start_time;
  clock_t stop_time;

  void start();
  
  void stop();

  Real elapsed() const;
  
  void print_elapsed(std::string preamble = "Elapsed time is ",
		     Real tare = 0.0) const;
};

#endif
