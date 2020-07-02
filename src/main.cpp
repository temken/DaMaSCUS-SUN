#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>
#include <mpi.h>

// Headers from libphysica
#include "Natural_Units.hpp"

// Headers from obscura
#include "DM_Distribution.hpp"
#include "DM_Particle_Standard.hpp"

#include "Data_Generation.hpp"
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

	//Initial terminal output
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
	MPI_Barrier(MPI_COMM_WORLD);
	////////////////////////////////////////////////////////////////////////

	Solar_Model SSM;

	obscura::DM_Particle_SI DM(0.1 * GeV);
	DM.Set_Sigma_Proton(1e-1 * pb);
	DM.Set_Sigma_Electron(0.0 * pb);
	DM.Print_Summary(mpi_rank);

	unsigned int sample_size		 = 10;
	double u_min					 = 0.0;
	unsigned int isoreflection_rings = 1;
	Simulation_Data data_set(sample_size, u_min, isoreflection_rings);
	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 1000, 50);
	data_set.Generate_Data(DM, SSM);
	data_set.Print_Summary(mpi_rank);

	obscura::Standard_Halo_Model SHM;
	Reflection_Spectrum spectrum(data_set, SSM, SHM, DM.mass);
	spectrum.Print_Summary(mpi_rank);

	////////////////////////////////////////////////////////////////////////
	//Final terminal output
	MPI_Barrier(MPI_COMM_WORLD);
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	if(mpi_rank == 0)
	{
		std::cout << "\n[Finished in " << std::round(1000. * durationTotal) / 1000. << "s";
		if(durationTotal > 60.0)
			std::cout << " (" << floor(durationTotal / 3600.0) << ":" << floor(fmod(durationTotal / 60.0, 60.0)) << ":" << floor(fmod(durationTotal, 60.0)) << ")]." << std::endl;
		else
			std::cout << "]" << std::endl;
	}
	MPI_Finalize();
	return 0;
}