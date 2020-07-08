#include "Parameter_Scan.hpp"

#include <algorithm>

// Headers from libphysica
#include "Utilities.hpp"

#include "Data_Generation.hpp"
#include "Reflection_Spectrum.hpp"

Parameter_Scan::Parameter_Scan(const std::vector<double> masses, const std::vector<double>& coupl, unsigned int samplesize)
: DM_masses(masses), couplings(coupl), sample_size(samplesize)
{
	std::sort(DM_masses.begin(), DM_masses.end());
	std::sort(couplings.begin(), couplings.end());
	p_value_grid = std::vector<std::vector<double>>(DM_masses.size(), std::vector<double>(couplings.size(), 1.0));
}

void Parameter_Scan::Perform_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, int mpi_rank)
{
	double mDM_original		 = DM.mass;
	double coupling_original = DM.Get_Interaction_Parameter("Electrons");

	obscura::Standard_Halo_Model SHM;
	unsigned int counter	  = 0;
	unsigned int mass_counter = 0;
	for(unsigned int i = 0; i < DM_masses.size(); i++)
	{
		mass_counter++;
		DM.Set_Mass(DM_masses[i]);

		unsigned int cs_counter = 0;
		for(unsigned int j = 0; j < couplings.size(); j++)
		{
			counter++;
			cs_counter++;
			if(mpi_rank == 0)
				std::cout << std::endl
						  << counter << ") Mass " << mass_counter << " / " << DM_masses.size() << "\t Cross section " << cs_counter << std::endl;

			DM.Set_Interaction_Parameter(couplings[couplings.size() - 1 - j], "Electrons");
			solar_model.Interpolate_Total_DM_Scattering_Rate(DM, 1000, 50);
			double u_min = detector.Minimum_DM_Speed(DM);
			Simulation_Data data_set(sample_size, u_min);
			data_set.Generate_Data(DM, solar_model);
			// data_set.Print_Summary(mpi_rank);
			Reflection_Spectrum spectrum(data_set, solar_model, SHM, DM.mass);
			spectrum.Print_Summary(mpi_rank);
			double p								  = detector.P_Value(DM, spectrum);
			p_value_grid[i][couplings.size() - 1 - j] = p;
			if(mpi_rank == 0)
				std::cout << "p-value = " << libphysica::Round(p) << std::endl;
			if(p > 0.25)
				break;
		}
	}
	DM.Set_Mass(mDM_original);
	DM.Set_Interaction_Parameter(coupling_original, "Electrons");
}

void Parameter_Scan::Export_P_Values(const std::string& path)
{
	std::vector<std::vector<double>> table;
	for(unsigned int i = 0; i < DM_masses.size(); i++)
		for(unsigned int j = 0; j < couplings.size(); j++)
			table.push_back({DM_masses[i], couplings[j], p_value_grid[i][j]});
	libphysica::Export_Table(path, table);
}