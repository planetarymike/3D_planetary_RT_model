//my_clock.h

#ifndef __my_clock_H
#define __my_clock_H

#include <ctime>
#include <iostream>

struct my_clock {
  clock_t start_time;
  clock_t stop_time;

  void start() {
    start_time=clock();
  }
  
  void stop() {
    stop_time = clock();
  }

  void print_elapsed() {
    double secs = (stop_time-start_time)*1.0/CLOCKS_PER_SEC;
    int hrs = secs/3600;
    int mins = (secs - hrs*3600)/60;
    secs = secs - hrs*3600 - mins*60;
    std::cout << "Elapsed time is " 
	      << hrs << " hours, " 
	      << mins << " minutes, and " 
	      << secs << " seconds.\n";
  }
};

#endif
