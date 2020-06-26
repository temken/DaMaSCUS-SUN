#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>
#include <random>

// Headers from libphysica
#include "Natural_Units.hpp"

// Headers from obscura
#include "Astronomy.hpp"
#include "DM_Distribution.hpp"
#include "DM_Particle_Standard.hpp"

#include "Data_Generation.hpp"
#include "Reflection_Spectrum.hpp"
#include "Simulation_Trajectory.hpp"
#include "Simulation_Utilities.hpp"
#include "Solar_Model.hpp"
#include "version.hpp"

using namespace libphysica::natural_units;

int main()
{
	//Initial terminal output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	std::cout << "[Started on " << ctime_start << "]" << std::endl;
	std::cout << PROJECT_NAME << "-" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
			  << DAMASCUS_SUN_LOGO << std::endl;
	////////////////////////////////////////////////////////////////////////

	std::random_device rd;
	std::mt19937 PRNG(rd());

	Solar_Model SSM;
	SSM.Print_Summary();

	// obscura::Standard_Halo_Model SHM;
	// SHM.Set_Observer_Velocity(libphysica::Vector({0, 0, 0}));

	obscura::DM_Particle_SI DM(0.25 * GeV);
	DM.Set_Sigma_Proton(1e-1 * pb);
	DM.Set_Sigma_Electron(1e-1 * pb);
	DM.Print_Summary();

	unsigned int sample_size		 = 1;
	double u_min					 = 0.0;
	unsigned int isoreflection_rings = 1;
	Simulation_Data data_set(sample_size, u_min, isoreflection_rings);
	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 1000, 50);
	data_set.Generate_Data(DM, SSM);
	data_set.Print_Summary();

	Reflection_Spectrum spectrum(data_set);
	spectrum.Print_Summary();

	////////////////////////////////////////////////////////////////////////
	//Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << std::round(1000. * durationTotal) / 1000. << "s";
	if(durationTotal > 60.0)
		std::cout << " (" << floor(durationTotal / 3600.0) << ":" << floor(fmod(durationTotal / 60.0, 60.0)) << ":" << floor(fmod(durationTotal, 60.0)) << ")]." << std::endl;
	else
		std::cout << "]" << std::endl;

	return 0;
}