#include <iostream>
#include <chrono>
#include <cmath>

#include "module.hpp"

// Headers from libphysica
#include "Natural_Units.hpp"

// Headers from obscura
#include "Astronomy.hpp"

int main()
{
	//Starting time
	auto time_start = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	std::cout 	<<"Started at " <<std::ctime(&time_start_t)<<std::endl;
	
	std::cout <<fib(10) <<std::endl;
	std::cout <<In_Units(1.0, meter/sec)<<std::endl;
	std::cout <<Fractional_Days_since_J2000(1,1,2001,12)<<std::endl;

	//Ending time and computing time
	auto time_end = std::chrono::system_clock::now();
	double durationTotal =1e-6*std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start ).count();
	std::cout 	<<"\n[Finished in "<< std::round(1000.*durationTotal)/1000.<<"s";
	if(durationTotal > 60.0)
		std::cout <<" ("<<floor(durationTotal/3600.0)<<":"<<floor(fmod(durationTotal/60.0,60.0))<<":"<<floor(fmod(durationTotal,60.0))<<")]."<<std::endl;
	else 
		std::cout <<"]"<<std::endl;
	
	return 0;
}