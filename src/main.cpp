#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>
#include <mpi.h>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "Data_Generation.hpp"
#include "Parameter_Scan.hpp"
#include "Reflection_Spectrum.hpp"
#include "Scattering_Rates.hpp"
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
	Solar_Model SSM(cfg.use_medium_effects, cfg.zeta);
	cfg.Print_Summary(mpi_rank);
	SSM.Print_Summary(mpi_rank);
	MPI_Barrier(MPI_COMM_WORLD);
	////////////////////////////////////////////////////////////////////////

	// Generate data for one parameter point specified in the configuration file.
	if(cfg.run_mode == "Parameter point")
	{
		double u_min = 0.0;
		// double u_min = cfg.DM_detector->Minimum_DM_Speed(*cfg.DM);
		Simulation_Data data_set(cfg.sample_size, u_min, cfg.isoreflection_rings);
		data_set.Configure(1.1 * rSun, 1, 1000);
		if(mpi_rank == 0)
			std::cout << "\nDM parameters:" << std::endl
					  << "\tm_DM [MeV]:\t" << libphysica::Round(In_Units(cfg.DM->mass, MeV)) << std::endl
					  << "\tsigma_p [cm2]:\t" << libphysica::Round(In_Units(cfg.DM->Get_Interaction_Parameter("Nuclei"), cm * cm)) << std::endl
					  << "\tsigma_e [cm2]:\t" << libphysica::Round(In_Units(cfg.DM->Get_Interaction_Parameter("Electrons"), cm * cm)) << std::endl
					  << std::endl;
		SSM.Interpolate_Total_DM_Scattering_Rate(*cfg.DM, cfg.interpolation_points, cfg.interpolation_points, mpi_rank);
		if(mpi_rank == 0)
			std::cout << "\nGenerating data...\t(with u_min = " << libphysica::Round(In_Units(u_min, km / sec)) << " km/s)" << std::endl;
		data_set.Generate_Data(*cfg.DM, SSM, *cfg.DM_distr);
		data_set.Print_Summary(mpi_rank);
		if(cfg.isoreflection_rings == 1)
		{
			Reflection_Spectrum spectrum(data_set, SSM, *cfg.DM_distr, cfg.DM->mass, 0);
			spectrum.Print_Summary(mpi_rank);

			// Export differential DM flux dPhi/dv to file
			std::function<double(double)> dPhi_dv = [&spectrum, &cfg](double v) {
				return spectrum.Differential_DM_Flux(v, cfg.DM->mass);
			};
			std::vector<double> speeds = libphysica::Linear_Space(spectrum.Minimum_DM_Speed(), spectrum.Maximum_DM_Speed(), 300);
			if(mpi_rank == 0)
				libphysica::Export_Function(cfg.results_path + "Differential_SRDM_Flux.txt", dPhi_dv, speeds, {km / sec, 1.0 / (km / sec) / cm / cm / sec});

			// Export recoil energy spectrum dR/dE to file
			std::function<double(double)> dR_dE = [&spectrum, &cfg](double E) {
				return cfg.DM_detector->dRdE(E, *cfg.DM, spectrum);
			};
			std::vector<double> energies = libphysica::Log_Space(0.1 * eV, cfg.DM_detector->Maximum_Energy_Deposit(*cfg.DM, spectrum), 300);
			if(mpi_rank == 0)
				libphysica::Export_Function(cfg.results_path + "Differential_Energy_Spectrum.txt", dR_dE, energies, {keV, 1.0 / keV / kg / year});

			// Compute p value for chosen experiment
			double p = cfg.DM_detector->P_Value(*cfg.DM, spectrum);
			libphysica::Print_Box("p = " + std::to_string(libphysica::Round(p)), 1, mpi_rank);
		}
		else
		{
			std::ofstream f;
			f.open(cfg.results_path + "/Detection_Rate.txt");
			std::vector<double> isoreflection_angles = Isoreflection_Ring_Angles(cfg.isoreflection_rings);
			if(mpi_rank == 0)
				std::cout << "Theta [deg]\t<u> [km/sec]\tDM flux [cm^-2 sec^-1]\tSignal rate [g^-1 day^-1]" << std::endl;

			for(unsigned int ring = 0; ring < cfg.isoreflection_rings; ring++)
			{
				Reflection_Spectrum spectrum(data_set, SSM, *cfg.DM_distr, cfg.DM->mass, ring);
				std::function<double(double)> func = [&spectrum, &cfg](double v) {
					return spectrum.Differential_DM_Flux(v, cfg.DM->mass);
				};
				std::vector<double> speeds = libphysica::Linear_Space(spectrum.Minimum_DM_Speed(), spectrum.Maximum_DM_Speed(), 300);
				if(mpi_rank == 0)
				{
					libphysica::Export_Function(cfg.results_path + "Differential_SRDM_Flux_" + std::to_string(ring) + ".txt", func, speeds, {km / sec, 1.0 / (km / sec) / cm / cm / sec});
					double total_rate = cfg.DM_detector->DM_Signal_Rate_Total(*cfg.DM, spectrum);
					f << isoreflection_angles[ring] << "\t" << spectrum.Average_Speed() / km * sec << "\t" << In_Units(spectrum.Total_DM_Flux(cfg.DM->mass), 1.0 / cm / cm / sec) << "\t" << In_Units(total_rate, 1.0 / gram / day) << std::endl;
					std::cout << libphysica::Round(isoreflection_angles[ring] / deg) << "\t\t" << libphysica::Round(In_Units(spectrum.Average_Speed(), km / sec)) << "\t\t" << libphysica::Round(In_Units(spectrum.Total_DM_Flux(cfg.DM->mass), 1.0 / cm / cm / sec)) << "\t\t\t" << libphysica::Round(In_Units(total_rate, 1.0 / gram / day)) << std::endl;
				}
			}
			f.close();
		}
	}
	// Perform a parameter scan to compute exclusion limits
	else if(cfg.run_mode == "Parameter scan")
	{
		if(mpi_rank == 0 && cfg.compute_halo_constraints)
		{
			std::cout << "\nCompute halo constraints for " << cfg.DM_detector->name << ":" << std::endl;
			double mDM_min								= cfg.DM_detector->Minimum_DM_Mass(*cfg.DM, *cfg.DM_distr);
			std::vector<double> DM_masses				= libphysica::Log_Space(mDM_min, GeV, 100);
			std::vector<std::vector<double>> halo_limit = cfg.DM_detector->Upper_Limit_Curve(*cfg.DM, *cfg.DM_distr, DM_masses, cfg.constraints_certainty);
			int CL										= std::round(100.0 * cfg.constraints_certainty);
			libphysica::Export_Table(TOP_LEVEL_DIR "results/" + cfg.ID + "/Halo_Limit_" + std::to_string(CL) + ".txt", halo_limit, {GeV, cm * cm});
		}
		Parameter_Scan scan(cfg);
		if(cfg.perform_full_scan)
			scan.Perform_Full_Scan(*cfg.DM, *cfg.DM_detector, SSM, *cfg.DM_distr, mpi_rank);
		else
			scan.Perform_STA_Scan(*cfg.DM, *cfg.DM_detector, SSM, *cfg.DM_distr, mpi_rank);
		scan.Export_Results(mpi_rank);
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
		std::random_device rd;
		std::mt19937 PRNG(rd());
		double r				 = 0.1 * rSun;
		double temperature		 = SSM.Temperature(r);
		int target_index		 = 4;
		obscura::Isotope nucleus = SSM.target_isotopes[target_index];
		std::cout << nucleus.name << std::endl;
		double nucleus_density						= SSM.Number_Density_Nucleus(target_index, r);
		double number_density_electrons				= SSM.Number_Density_Electron(r);
		auto nuclei									= SSM.target_isotopes;
		std::vector<double> number_densities_nuclei = SSM.Number_Densities_Nuclei(r);
		double vDM									= 5000 * km / sec;
		double qMin									= SSM.zeta * cfg.DM->mass * vDM;
		double qMax_Nucleus							= Maximum_Momentum_Transfer(cfg.DM->mass, temperature, nucleus.mass, vDM, 5);
		double qMax_Electron						= Maximum_Momentum_Transfer(cfg.DM->mass, temperature, mElectron, vDM, 5);
		auto qList_Nucleus							= libphysica::Linear_Space(qMin, qMax_Nucleus, 1000);
		auto qList_Electron							= libphysica::Linear_Space(qMin, qMax_Electron, 1000);

		std::ofstream f, g;
		// // 0. Differential scattering rate
		// int nPoints	 = 200;
		// auto qList	 = libphysica::Linear_Space(qMin, qMax, nPoints);
		// auto cosList = libphysica::Linear_Space(-1, 0, nPoints);
		// f.open("Differential_Scattering_Rate.txt");
		// for(auto q : qList)
		// 	for(auto cos : cosList)
		// 	{
		// 		double differential_scattering_rate = Differential_Scattering_Rate_Nucleus(q, cos, *cfg.DM, vDM, nucleus, nucleus_density, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects);
		// 		f << q / qMax << "\t" << cos << "\t" << differential_scattering_rate << std::endl;
		// 	}
		// f.close();

		// 1. PDF of cos(theta)
		std::cout << "1. PDF of cos(theta)" << std::endl;
		f.open("PDF_Cos_Theta_Electron.txt");
		g.open("PDF_Cos_Theta_Nucleus.txt");
		for(auto& x : libphysica::Linear_Space(-1.0, 1.0, 1000))
		{
			f << x << "\t" << PDF_Cos_Theta_Electron(x, *cfg.DM, vDM, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects, qMin) << std::endl;
			g << x << "\t" << PDF_Cos_Theta_Nucleus(x, *cfg.DM, vDM, nucleus, nucleus_density, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects, qMin) << std::endl;
		}
		f.close();
		g.close();

		// 2. Sample cos(theta)
		std::cout << "2. Sample cos(theta)" << std::endl;
		f.open("Sample_Cos_Theta_Electron.txt");
		g.open("Sample_Cos_Theta_Nucleus.txt");
		for(int i = 0; i < cfg.sample_size; i++)
		{
			f << Sample_Cos_Theta_Electron(PRNG, *cfg.DM, vDM, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects, qMin) << std::endl;
			g << Sample_Cos_Theta_Nucleus(PRNG, *cfg.DM, vDM, nucleus, nucleus_density, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects, qMin) << std::endl;
		}
		f.close();
		g.close();

		// 3. PDF of q
		std::cout << "3. PDF of q" << std::endl;
		f.open("PDF_q_Electron.txt");
		g.open("PDF_q_Nucleus.txt");
		double cos_theta_e = Sample_Cos_Theta_Electron(PRNG, *cfg.DM, vDM, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects, qMin);
		double cos_theta_n = Sample_Cos_Theta_Nucleus(PRNG, *cfg.DM, vDM, nucleus, nucleus_density, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects, qMin);
		for(unsigned int i = 0; i < qList_Nucleus.size(); i++)
		{
			f << qList_Electron[i] << "\t" << PDF_q_Electron(qList_Electron[i], cos_theta_e, *cfg.DM, vDM, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects, qMin) << std::endl;
			g << qList_Nucleus[i] << "\t" << PDF_q_Nucleus(qList_Nucleus[i], cos_theta_n, *cfg.DM, vDM, nucleus, nucleus_density, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects, qMin) << std::endl;
		}
		f.close();
		g.close();

		// 4. CDF of q
		std::cout << "4. CDF of q" << std::endl;
		f.open("CDF_q_Electron.txt");
		g.open("CDF_q_Nucleus.txt");
		for(unsigned int i = 0; i < qList_Nucleus.size(); i++)

		{
			f << qList_Electron[i] << "\t" << CDF_q_Electron(qList_Electron[i], cos_theta_e, *cfg.DM, vDM, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects, qMin) << std::endl;
			g << qList_Nucleus[i] << "\t" << CDF_q_Nucleus(qList_Nucleus[i], cos_theta_n, *cfg.DM, vDM, nucleus, nucleus_density, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects, qMin) << std::endl;
		}
		f.close();
		g.close();

		// 5. Sample q
		std::cout << "5. Sample q" << std::endl;
		f.open("Sample_q_Electron.txt");
		g.open("Sample_q_Nucleus.txt");
		for(int i = 0; i < cfg.sample_size; i++)
		{
			f << Sample_q_Electron(PRNG, cos_theta_e, *cfg.DM, vDM, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects, qMin) << std::endl;
			g << Sample_q_Nucleus(PRNG, cos_theta_n, *cfg.DM, vDM, nucleus, nucleus_density, temperature, number_density_electrons, nuclei, number_densities_nuclei, SSM.use_medium_effects, qMin) << std::endl;
		}
		f.close();
		g.close();
	}

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	MPI_Barrier(MPI_COMM_WORLD);
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	if(mpi_rank == 0)
		std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	MPI_Finalize();
	return 0;
}