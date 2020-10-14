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
	Solar_Model SSM;
	cfg.Print_Summary(mpi_rank);
	MPI_Barrier(MPI_COMM_WORLD);
	////////////////////////////////////////////////////////////////////////

	// Generate data for one parameter point specified in the configuration file.
	if(cfg.run_mode == "Parameter point")
	{
		SSM.Interpolate_Total_DM_Scattering_Rate(*cfg.DM, 1000, 1000);
		double u_min = 0.0;
		// double u_min = cfg.DM_detector->Minimum_DM_Speed(*cfg.DM);
		Simulation_Data data_set(cfg.sample_size, u_min, cfg.isoreflection_rings);
		data_set.Configure(1.1 * rSun, 1, 1000);
		if(mpi_rank == 0)
			std::cout << "Generate data..." << std::endl
					  << "\tm_DM [MeV]:\t" << libphysica::Round(In_Units(cfg.DM->mass, MeV)) << "\t\t"
					  << "sigma_p [cm2]:\t" << libphysica::Round(In_Units(cfg.DM->Get_Interaction_Parameter("Nuclei"), cm * cm)) << std::endl
					  << "\tu_min [km/sec]:\t" << libphysica::Round(In_Units(u_min, km / sec)) << "\t\t"
					  << "sigma_e [cm2]:\t" << libphysica::Round(In_Units(cfg.DM->Get_Interaction_Parameter("Electrons"), cm * cm)) << std::endl
					  << std::endl;
		data_set.Generate_Data(*cfg.DM, SSM, *cfg.DM_distr);
		data_set.Print_Summary(mpi_rank);
		Reflection_Spectrum spectrum(data_set, SSM, *cfg.DM_distr, cfg.DM->mass, 0);
		spectrum.Print_Summary(mpi_rank);
		double p = cfg.DM_detector->P_Value(*cfg.DM, spectrum);
		libphysica::Print_Box("p = " + std::to_string(libphysica::Round(p)), 1, mpi_rank);
	}
	//Perform a parameter scan to compute exclusion limits
	else if(cfg.run_mode == "Parameter scan")
	{
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
		scan.Perform_STA_Scan(*cfg.DM, *cfg.DM_detector, SSM, *cfg.DM_distr, cfg.ID, mpi_rank);
		scan.Export_Results(cfg.ID, mpi_rank);
		if(mpi_rank == 0)
		{
			int CL = std::round(100.0 * cfg.constraints_certainty);
			std::cout << "\nFinal reflection constraints (" << CL << "% CL)" << std::endl;
			scan.Print_Grid(mpi_rank);
		}
	}
	// Run some custom code
	else
	{
		// Test isotropy
		if(mpi_rank == 0)
			std::cout << "Test isotropy of reflection flux." << std::endl;
		std::ofstream f;
		f.open(TOP_LEVEL_DIR "results/" + cfg.ID + "/Anisotropy.txt");
		std::vector<double> angles = Isoreflection_Ring_Angles(cfg.isoreflection_rings);
		SSM.Interpolate_Total_DM_Scattering_Rate(*cfg.DM, 1000, 1000);
		double u_min = 0.0;
		// std::cout << cfg.DM_detector->Minimum_DM_Speed(*cfg.DM) / km * sec << std::endl;
		Simulation_Data data_set(cfg.sample_size, u_min, cfg.isoreflection_rings);
		if(mpi_rank == 0)
			std::cout << "Generate data..." << std::endl
					  << "\tm_DM [MeV]:\t" << libphysica::Round(In_Units(cfg.DM->mass, MeV)) << "\t\t"
					  << "sigma_p [cm2]:\t" << libphysica::Round(In_Units(cfg.DM->Get_Interaction_Parameter("Nuclei"), cm * cm)) << std::endl
					  << "\tu_min [km/sec]:\t" << libphysica::Round(In_Units(u_min, km / sec)) << "\t\t"
					  << "sigma_e [cm2]:\t" << libphysica::Round(In_Units(cfg.DM->Get_Interaction_Parameter("Electrons"), cm * cm)) << std::endl
					  << std::endl;

		data_set.Generate_Data(*cfg.DM, SSM, *cfg.DM_distr);
		data_set.Print_Summary(mpi_rank);
		if(mpi_rank == 0)
		{
			for(unsigned int ring = 0; ring < cfg.isoreflection_rings; ring++)
			{
				Reflection_Spectrum spectrum(data_set, SSM, *cfg.DM_distr, cfg.DM->mass, ring);
				double total_rate = cfg.DM_detector->DM_Signals_Total(*cfg.DM, spectrum) / (0.1 * kg * year);
				f << angles[ring] << "\t" << spectrum.Average_Speed() / km * sec << "\t" << In_Units(spectrum.Total_DM_Flux(cfg.DM->mass), 1.0 / cm / cm / sec) << "\t" << In_Units(total_rate, 1.0 / gram / day) << std::endl;			  //<< "\t" << In_Units(total_rate_Xe, 1.0 / gram / day) << std::endl;
				std::cout << angles[ring] << "\t" << spectrum.Average_Speed() / km * sec << "\t" << In_Units(spectrum.Total_DM_Flux(cfg.DM->mass), 1.0 / cm / cm / sec) << "\t" << In_Units(total_rate, 1.0 / gram / day) << std::endl;	  //<< "\t" << In_Units(total_rate_Xe, 1.0 / gram / day) << std::endl;
			}
		}
		f.close();
	}

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