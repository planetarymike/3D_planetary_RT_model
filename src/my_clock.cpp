//my_clock.cpp

#include "my_clock.hpp"

void my_clock::start() {
  start_time=clock();
}
  
void my_clock::stop() {
  stop_time = clock();
}

Real my_clock::elapsed() const {
  return (stop_time-start_time)*1.0/CLOCKS_PER_SEC;
}
  
void my_clock::print_elapsed(std::string preamble,
			     Real tare) const {
  Real secs = elapsed() - tare;

  std::cout << preamble;
  if (secs < 0.001) {
    std::string mu = "\u03BC";
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
