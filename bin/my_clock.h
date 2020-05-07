//my_clock.h

#ifndef __my_clock_H
#define __my_clock_H

#include <ctime>
#include <iostream>
#include <string>

using std::string;

struct my_clock {
  clock_t start_time;
  clock_t stop_time;

  void start() {
    start_time=clock();
  }
  
  void stop() {
    stop_time = clock();
  }

  double elapsed() {
    return (stop_time-start_time)*1.0/CLOCKS_PER_SEC;
  }
  
  void print_elapsed(string preamble = "Elapsed time is ", double tare = 0.0) {
    double secs = elapsed() - tare;

    std::cout << preamble;
    if (secs < 0.001) {
      string mu = "\u03BC";
      std::cout << (int) (secs*1000000) << " " << mu << "s .\n";
    } else if (secs < 1) {
      std::cout << (int) (secs*1000) << " ms .\n";
    } else if (secs < 60) {
      std::cout << secs << " s .\n";
    } else if (secs < 3600) {
      int mins = secs/60;
      secs = secs - mins*60;
      std::cout << mins << " minutes, and " 
		<< secs << " seconds.\n";
      
    } else {
      int hrs = secs/3600;
      int mins = (secs - hrs*3600)/60;
      secs = secs - hrs*3600 - mins*60;
      std::cout << hrs << " hours, " 
		<< mins << " minutes, and " 
		<< secs << " seconds.\n";
    }
  }
};

#endif
