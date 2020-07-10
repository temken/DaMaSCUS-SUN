#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>
#include <mpi.h>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Utilities.hpp"

// Headers from obscura
#include "Astronomy.hpp"
#include "Configuration.hpp"
#include "DM_Distribution.hpp"
#include "DM_Particle_Standard.hpp"
#include "Experiments.hpp"

#include "Data_Generation.hpp"
#include "Parameter_Scan.hpp"
#include "Reflection_Spectrum.hpp"
#include "Solar_Model.hpp"
#include "version.hpp"

using namespace libphysica::natural_units;

int main()
{
	MPI_Init(NULL, NULL);
	int mpi_processes, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	// Initial terminal output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	if(mpi_rank == 0)
		std::cout << "[Started on " << ctime_start << "]" << std::endl
				  << PROJECT_NAME << "-" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
				  << DAMASCUS_SUN_LOGO
				  << std::endl
				  << "MPI processes:\t" << mpi_processes << std::endl;

	// Configuration parameters
	obscura::Configuration cfg(PROJECT_DIR "bin/config.cfg", mpi_rank);
	Solar_Model SSM;
	// cfg.Print_Summary(mpi_rank);
	MPI_Barrier(MPI_COMM_WORLD);
	////////////////////////////////////////////////////////////////////////

	// ////////////////////////////////////////////////////////////////////////
	// // Parameter Scan
	// std::vector<double> DM_masses	   = libphysica::Log_Space(1.0 * keV, 5.0 * MeV, 20);
	// std::vector<double> cross_sections = libphysica::Log_Space(1.0e-35 * cm * cm, 1.0e-32 * cm * cm, 15);
	// unsigned int sample_size		   = 10;
	// Parameter_Scan scan(DM_masses, cross_sections, sample_size);
	// scan.Perform_Scan(*cfg.DM, *cfg.DM_detector, SSM, *cfg.DM_distr, mpi_rank);
	// scan.Export_P_Values(TOP_LEVEL_DIR "results/" + cfg.ID + "/p_values_3.txt");
	// // scan.Import_P_Values(TOP_LEVEL_DIR "results/" + cfg.ID + "/p_values.txt");
	// scan.Export_Limits(TOP_LEVEL_DIR "results/" + cfg.ID + "/", mpi_rank);
	// ////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////
	// Generate data
	SSM.Interpolate_Total_DM_Scattering_Rate(*cfg.DM, 1000, 50);
	unsigned int sample_size = 100;
	double u_min			 = 0.0;
	Simulation_Data data_set(sample_size, u_min);
	// data_set.Configure(1.1 * rSun, 0, 100);
	data_set.Generate_Data(*cfg.DM, SSM, *cfg.DM_distr);
	data_set.Print_Summary(mpi_rank);
	Reflection_Spectrum spectrum(data_set, SSM, *cfg.DM_distr, cfg.DM->mass);
	spectrum.Print_Summary(mpi_rank);
	////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////
	//Final terminal output
	MPI_Barrier(MPI_COMM_WORLD);
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	if(mpi_rank == 0)
		std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]" << std::endl;
	MPI_Finalize();
	return 0;
}