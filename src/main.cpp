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

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

int main(int argc, char* argv[])
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
	Configuration cfg(argv[1], mpi_rank);
	// Configuration cfg(PROJECT_DIR "bin/config.cfg", mpi_rank);
	Solar_Model SSM;
	cfg.Print_Summary(mpi_rank);
	MPI_Barrier(MPI_COMM_WORLD);
	////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////
	// Parameter Scan
	if(mpi_rank == 0 && cfg.compute_halo_constraints)
	{
		std::cout << "Compute halo constraints for " << cfg.DM_detector->name << ":" << std::endl;
		double mDM_min								= cfg.DM_detector->Minimum_DM_Mass(*cfg.DM, *cfg.DM_distr);
		std::vector<double> DM_masses				= libphysica::Log_Space(mDM_min, GeV, 100);
		std::vector<std::vector<double>> halo_limit = cfg.DM_detector->Upper_Limit_Curve(*cfg.DM, *cfg.DM_distr, DM_masses, cfg.constraints_certainty);
		int CL										= std::round(100.0 * cfg.constraints_certainty);
		libphysica::Export_Table(TOP_LEVEL_DIR "results/" + cfg.ID + "/Halo_Limit_" + std::to_string(CL) + ".txt", halo_limit, {GeV, cm * cm});
	}
	Parameter_Scan scan(cfg);
	// scan.Perform_Full_Scan(*cfg.DM, *cfg.DM_detector, SSM, *cfg.DM_distr, mpi_rank);
	// scan.Import_P_Values(cfg.ID);
	scan.Perform_STA_Scan(*cfg.DM, *cfg.DM_detector, SSM, *cfg.DM_distr, mpi_rank);
	scan.Export_Results(cfg.ID, mpi_rank);
	if(mpi_rank == 0)
	{
		int CL = std::round(100.0 * cfg.constraints_certainty);
		std::cout << "\nFinal reflection constraints (" << CL << "% CL)" << std::endl;
		scan.Print_Grid(mpi_rank);
	}
	// ////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////
	// Generate data
	// SSM.Interpolate_Total_DM_Scattering_Rate(*cfg.DM, 1000, 50);
	// unsigned int sample_size = 1000;
	// double u_min			 = 0.0;
	// // double u_min = cfg.DM_detector->Minimum_DM_Speed(*cfg.DM);
	// Simulation_Data data_set(sample_size, u_min);
	// // data_set.Configure(1.1 * rSun, 1, 1);
	// data_set.Generate_Data(*cfg.DM, SSM, *cfg.DM_distr);
	// data_set.Print_Summary(mpi_rank);
	// Reflection_Spectrum spectrum(data_set, SSM, *cfg.DM_distr, cfg.DM->mass);
	// spectrum.Print_Summary(mpi_rank);
	// double p = cfg.DM_detector->P_Value(*cfg.DM, spectrum);
	// if(mpi_rank == 0)
	// 	std::cout << "p-value = " << libphysica::Round(p) << std::endl;
	////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////
	//Final terminal output
	MPI_Barrier(MPI_COMM_WORLD);
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	if(mpi_rank == 0)
		std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	MPI_Finalize();
	return 0;
}