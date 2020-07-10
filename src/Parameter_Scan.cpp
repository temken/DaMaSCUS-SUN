#include "Parameter_Scan.hpp"

#include <algorithm>
#include <set>

// Headers from libphysica
#include "Utilities.hpp"

#include "Data_Generation.hpp"
#include "Reflection_Spectrum.hpp"

Parameter_Scan::Parameter_Scan(const std::vector<double> masses, const std::vector<double>& coupl, unsigned int samplesize)
: DM_masses(masses), couplings(coupl), sample_size(samplesize)
{
	std::sort(DM_masses.begin(), DM_masses.end());
	std::sort(couplings.begin(), couplings.end());
	p_value_grid = std::vector<std::vector<double>>(couplings.size(), std::vector<double>(DM_masses.size(), 1.0));
}

void Parameter_Scan::Perform_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, int mpi_rank)
{
	double mDM_original		 = DM.mass;
	double coupling_original = DM.Get_Interaction_Parameter(detector.Target_Particles());

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

			DM.Set_Interaction_Parameter(couplings[couplings.size() - 1 - j], detector.Target_Particles());
			solar_model.Interpolate_Total_DM_Scattering_Rate(DM, 1000, 50);
			double u_min = detector.Minimum_DM_Speed(DM);
			Simulation_Data data_set(sample_size, u_min);

			data_set.Generate_Data(DM, solar_model);
			Reflection_Spectrum spectrum(data_set, solar_model, SHM, DM.mass);

			double p								  = detector.P_Value(DM, spectrum);
			p_value_grid[couplings.size() - 1 - j][i] = p;
			if(mpi_rank == 0)
				std::cout << "p-value = " << libphysica::Round(p) << std::endl;
			if(p > 0.15)
				break;
		}
	}
	DM.Set_Mass(mDM_original);
	DM.Set_Interaction_Parameter(coupling_original, detector.Target_Particles());
}

std::vector<std::vector<double>> Parameter_Scan::Limit_Curve(double certainty_level)
{
	std::vector<std::vector<double>> limit;
	for(unsigned int i = 0; i < DM_masses.size(); i++)
	{
		if(p_value_grid.back()[i] < (1.0 - certainty_level))
		{
			std::vector<double> interpolation_list;
			for(unsigned int j = 0; j < couplings.size(); j++)
				interpolation_list.push_back(p_value_grid[j][i] - (1.0 - certainty_level));
			libphysica::Interpolation interpolation(couplings, interpolation_list);
			double coupling_limit = libphysica::Find_Root(interpolation, couplings[0], couplings.back(), 1.0e-6);
			limit.push_back({DM_masses[i], coupling_limit});
		}
	}
	return limit;
}

void Parameter_Scan::Import_P_Values(const std::string& file_path)
{
	std::vector<std::vector<double>> table = libphysica::Import_Table(file_path);
	// 1. Determine grid dimensions
	std::vector<double> all_masses;
	for(unsigned int i = 0; i < table.size(); i++)
		all_masses.push_back(table[i][0]);
	std::set<double> mass_set(all_masses.begin(), all_masses.end());
	unsigned int number_of_masses	 = DM_masses.size();
	unsigned int number_of_couplings = table.size() / number_of_masses;

	// 2. Assign masses and couplings
	DM_masses.assign(mass_set.begin(), mass_set.end());
	couplings.resize(number_of_couplings);
	for(unsigned int i = 0; i < number_of_couplings; i++)
		couplings[i] = table[i][1];
	//3. Assign the p values to the table
	unsigned int k = 0;
	p_value_grid   = std::vector<std::vector<double>>(number_of_couplings, std::vector<double>(number_of_masses, 0.0));
	for(unsigned int i = 0; i < number_of_masses; i++)
		for(unsigned int j = 0; j < number_of_couplings; j++)
			p_value_grid[j][i] = table[k++][2];
}

void Parameter_Scan::Export_P_Values(const std::string& folder_path, int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::vector<std::vector<double>> table;
		for(unsigned int i = 0; i < DM_masses.size(); i++)
			for(unsigned int j = 0; j < couplings.size(); j++)
				table.push_back({DM_masses[i], couplings[j], p_value_grid[j][i]});
		libphysica::Export_Table(folder_path + "p_values.txt", table);
	}
}

void Parameter_Scan::Export_Limits(const std::string& folder_path, int mpi_rank, std::vector<double> certainty_levels)
{
	for(auto certainty_level : certainty_levels)
	{
		std::vector<std::vector<double>> limit = Limit_Curve(certainty_level);
		int CL								   = std::round(100 * certainty_level);
		std::string filename				   = "Limit_" + std::to_string(CL) + ".txt";
		if(mpi_rank == 0)
			libphysica::Export_Table(folder_path + filename, limit);
	}
}