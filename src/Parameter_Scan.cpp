#include "Parameter_Scan.hpp"

#include <algorithm>
#include <libconfig.h++>
#include <set>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Utilities.hpp"

#include "Data_Generation.hpp"
#include "Reflection_Spectrum.hpp"

namespace DaMaSCUS_SUN
{

using namespace libconfig;
using namespace libphysica::natural_units;

Configuration::Configuration(std::string cfg_filename, int MPI_rank)
: obscura::Configuration(cfg_filename, MPI_rank)
{
	Import_Parameter_Scan_Parameter();
}

void Configuration::Import_Parameter_Scan_Parameter()
{
	try
	{
		sample_size = config.lookup("sample_size");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'sample_size' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		cross_section_min = config.lookup("cross_section_min");
		cross_section_min *= cm * cm;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'cross_section_min' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		cross_section_max = config.lookup("cross_section_max");
		cross_section_max *= cm * cm;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'cross_section_max' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		cross_sections = config.lookup("cross_sections");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'cross_sections' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		compute_halo_constraints = config.lookup("compute_halo_constraints");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'compute_halo_constraints' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void Configuration::Print_Summary(int mpi_rank)
{
	if(mpi_rank == 0)
	{
		Print_Summary_Base(mpi_rank);
		std::cout << "DaMaSCUS-SUN parameters" << std::endl
				  << "\tSample size:\t\t\t" << sample_size << std::endl
				  << "\tCross section (min) [cm^2]:\t" << libphysica::Round(In_Units(cross_section_min, cm * cm)) << std::endl
				  << "\tCross section (max) [cm^2]:\t" << libphysica::Round(In_Units(cross_section_max, cm * cm)) << std::endl
				  << "\tCross section steps:\t\t" << cross_sections << std::endl
				  << SEPARATOR << std::endl;
	}
}

Parameter_Scan::Parameter_Scan(Configuration& config)
: Parameter_Scan(libphysica::Log_Space(config.constraints_mass_min, config.constraints_mass_max, config.constraints_masses), libphysica::Log_Space(config.cross_section_min, config.cross_section_max, config.cross_sections), config.sample_size)
{
}

Parameter_Scan::Parameter_Scan(const std::vector<double>& masses, const std::vector<double>& coupl, unsigned int samplesize)
: DM_masses(masses), couplings(coupl), sample_size(samplesize)
{
	std::sort(DM_masses.begin(), DM_masses.end());
	std::sort(couplings.begin(), couplings.end());
	p_value_grid = std::vector<std::vector<double>>(couplings.size(), std::vector<double>(DM_masses.size(), -1.0));
}

void Parameter_Scan::Go_Forward(int& row, int& col, std::string& direction)
{
	if(direction == "N")

	{
		if(row == couplings.size() - 1)
			Go_Right(row, col, direction);
		else
			row++;
	}
	else if(direction == "E")
	{
		if(col == DM_masses.size() - 1)
			Go_Right(row, col, direction);
		else
			col++;
	}
	else if(direction == "S")
	{
		if(row == 0)
			Go_Right(row, col, direction);
		else
			row--;
	}
	else if(direction == "W")
	{
		if(col == 0)
			Go_Right(row, col, direction);
		else
			col--;
	}
}

void Parameter_Scan::Go_Left(int& row, int& col, std::string& direction)
{
	if(direction == "N")
	{
		if(col == 0)
			Go_Forward(row, col, direction);
		else
		{
			col--;
			direction = "W";
		}
	}
	else if(direction == "E")
	{
		if(row == couplings.size() - 1)
			Go_Forward(row, col, direction);
		else
		{
			row++;
			direction = "N";
		}
	}
	else if(direction == "S")
	{
		if(col == DM_masses.size() - 1)
			Go_Forward(row, col, direction);
		else
		{
			col++;
			direction = "E";
		}
	}
	else if(direction == "W")
	{
		if(row == 0)
			Go_Forward(row, col, direction);
		else
		{
			row--;
			direction = "S";
		}
	}
}

void Parameter_Scan::Go_Right(int& row, int& col, std::string& direction)
{
	if(direction == "N")
	{
		if(col == DM_masses.size() - 1)
			Go_Forward(row, col, direction);
		else
		{
			col++;
			direction = "E";
		}
	}
	else if(direction == "E")
	{
		if(row == 0)
			Go_Forward(row, col, direction);
		else
		{
			row--;
			direction = "S";
		}
	}
	else if(direction == "S")
	{
		if(col == 0)
			Go_Forward(row, col, direction);
		else
		{
			col--;
			direction = "W";
		}
	}
	else if(direction == "W")
	{
		if(row == couplings.size() - 1)
			Go_Forward(row, col, direction);
		else
		{
			row++;
			direction = "N";
		}
	}
}

void Parameter_Scan::Fill_STA_Gaps(double certainty_level)
{
	for(unsigned int row = 0; row < couplings.size(); row++)
	{
		bool excluded = false;
		for(unsigned int col = 0; col < DM_masses.size(); col++)
		{
			double p = p_value_grid[row][col];
			if(p < 0)
				p_value_grid[row][col] = excluded ? 0.0 : 1.0;
			else if(p < 1.0 - certainty_level)
				excluded = true;
			else
				excluded = false;
		}
	}
}

double Parameter_Scan::Compute_p_Value(int row, int col, obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank)
{
	if(p_value_grid[row][col] > 0)
		return p_value_grid[row][col];
	else
	{
		double mDM_original		 = DM.mass;
		double coupling_original = DM.Get_Interaction_Parameter(detector.Target_Particles());

		DM.Set_Interaction_Parameter(couplings[row], detector.Target_Particles());
		DM.Set_Mass(DM_masses[col]);
		double u_min = detector.Minimum_DM_Speed(DM);

		if(mpi_rank == 0)
			std::cout << std::endl
					  << ++counter << ")\t"
					  << "m_DM [MeV]:\t" << libphysica::Round(In_Units(DM.mass, MeV)) << "\t\t"
					  << "u_min [km/sec]:\t" << libphysica::Round(In_Units(u_min, km / sec)) << std::endl
					  << "\tsigma_p [cm2]:\t" << libphysica::Round(In_Units(DM.Get_Interaction_Parameter("Nuclei"), cm * cm)) << "\t\t"
					  << "sigma_e [cm2]:\t" << libphysica::Round(In_Units(DM.Get_Interaction_Parameter("Electrons"), cm * cm)) << std::endl
					  << std::endl;
		Print_Grid(mpi_rank, row, col);

		solar_model.Interpolate_Total_DM_Scattering_Rate(DM, 1000, 50);
		Simulation_Data data_set(sample_size, u_min);
		data_set.Generate_Data(DM, solar_model, halo_model);
		Reflection_Spectrum spectrum(data_set, solar_model, halo_model, DM.mass);
		double p = detector.P_Value(DM, spectrum);

		DM.Set_Mass(mDM_original);
		DM.Set_Interaction_Parameter(coupling_original, detector.Target_Particles());

		p_value_grid[row][col] = p;
		if(mpi_rank == 0)
		{
			std::cout << std::endl
					  << std::endl;
			libphysica::Print_Box("p = " + std::to_string(libphysica::Round(p)), 1);
		}
		return ((p < 1.0e-100) ? 0.0 : p);
	}
}

void Parameter_Scan::Perform_Full_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank)
{
	int last_excluded_mass_index = DM_masses.size();
	for(unsigned int i = 0; i < couplings.size(); i++)
	{
		int index_coupling = couplings.size() - 1 - i;
		bool row_exclusion = false;
		for(unsigned int j = 0; j < DM_masses.size(); j++)
		{
			int index_mass = DM_masses.size() - 1 - j;
			double p	   = Compute_p_Value(index_coupling, index_mass, DM, detector, solar_model, halo_model, mpi_rank);

			p_value_grid[index_coupling][index_mass] = p;

			if(mpi_rank == 0)
			{
				std::cout << std::endl
						  << std::endl;
				libphysica::Print_Box("p = " + std::to_string(libphysica::Round(p)), 1);
			}

			if(p < 0.1)
			{
				row_exclusion			 = true;
				last_excluded_mass_index = j;
			}
			else if(row_exclusion || j > last_excluded_mass_index + 1)
			{
				for(unsigned int k = 0; k < index_mass; k++)
					p_value_grid[index_coupling][k] = 1.0;
				break;
			}
		}
		if(!row_exclusion)
		{
			for(unsigned int k = 0; k < index_coupling; k++)
				for(unsigned int j = 0; j < DM_masses.size(); j++)
					p_value_grid[k][j] = 1.0;
			break;
		}
	}
}

void Parameter_Scan::Perform_STA_Scan(double certainty_level, obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank)
{
	std::string direction = "S";
	int row				  = couplings.size() - 1;
	int col				  = DM_masses.size() - 1;
	int row_start		  = row;
	int col_start		  = col;
	while(true)
	{
		double p = Compute_p_Value(row, col, DM, detector, solar_model, halo_model, mpi_rank);

		if(p < 1.0 - certainty_level)
			Go_Left(row, col, direction);
		else
		{
			if(col == 0 && direction == "S")
			{
				direction = "N";
				Go_Forward(row, col, direction);
			}
			else
				Go_Right(row, col, direction);
		}
		if(row == row_start && col == col_start)
			break;
	}
	Fill_STA_Gaps(certainty_level);
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
			double coupling_limit = libphysica::Find_Root(interpolation, couplings[0], couplings.back(), 0.01 * couplings[0]);
			limit.push_back({DM_masses[i], coupling_limit});
		}
	}
	return limit;
}

void Parameter_Scan::Import_P_Values(const std::string& ID)
{
	std::string filepath				   = TOP_LEVEL_DIR "results/" + ID + "/p_values.txt";
	std::vector<std::vector<double>> table = libphysica::Import_Table(filepath, {GeV, cm * cm, 1.0});
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

void Parameter_Scan::Export_P_Values(const std::string& ID, int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::vector<std::vector<double>> table;
		for(unsigned int i = 0; i < DM_masses.size(); i++)
			for(unsigned int j = 0; j < couplings.size(); j++)
				table.push_back({DM_masses[i], couplings[j], p_value_grid[j][i]});
		libphysica::Export_Table(TOP_LEVEL_DIR "results/" + ID + "/p_values.txt", table, {GeV, cm * cm, 1.0});
		libphysica::Export_Table(TOP_LEVEL_DIR "results/" + ID + "/p_grid.txt", p_value_grid);
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
			libphysica::Export_Table(folder_path + filename, limit, {GeV, cm * cm});
	}
}

void Parameter_Scan::Print_Grid(int mpi_rank, int marker_row, int marker_col)
{
	if(mpi_rank == 0)
	{
		for(int row = 0; row < couplings.size(); row++)
		{
			std::cout << "\t";
			for(int col = 0; col < DM_masses.size(); col++)
			{
				double p = p_value_grid[couplings.size() - 1 - row][col];
				if(row == couplings.size() - 1 - marker_row && col == marker_col)
					std::cout << "¤";
				else if(p < 0.0)
					std::cout << "·";
				else if(p < 0.1)
					std::cout << "█";
				else if(p > 0.1)
					std::cout << "░";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

Solar_Reflection_Limit::Solar_Reflection_Limit(unsigned int Nsample, double mMin, double mMax, unsigned int Nmass, double c_min, double c_max, double CL)
: sample_size(Nsample), coupling_min(c_min), coupling_max(c_max), certainty_level(CL)
{
	masses = libphysica::Log_Space(mMin, mMax, Nmass);
}

double Solar_Reflection_Limit::Upper_Limit(double mass, obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank)
{
	double mDM_original		 = DM.mass;
	double coupling_original = DM.Get_Interaction_Parameter(detector.Target_Particles());
	DM.Set_Mass(mass);
	std::function<double(double)> func = [this, &DM, &detector, &solar_model, &halo_model, mpi_rank](double log_coupling) {
		DM.Set_Interaction_Parameter(exp(log_coupling), detector.Target_Particles());
		solar_model.Interpolate_Total_DM_Scattering_Rate(DM, 1000, 50);
		double u_min = detector.Minimum_DM_Speed(DM);
		Simulation_Data data_set(sample_size, u_min);
		data_set.Generate_Data(DM, solar_model, halo_model);
		Reflection_Spectrum spectrum(data_set, solar_model, halo_model, DM.mass);
		double p = detector.P_Value(DM, spectrum);
		if(mpi_rank == 0)
			std::cout << "p = " << libphysica::Round(p) << std::endl;
		return p - (1.0 - certainty_level);
	};
	double log_coupling_min = log(coupling_min);
	double log_coupling_max = log(coupling_max);
	double log_limit		= libphysica::Find_Root(func, log_coupling_min, log_coupling_max, 1.0e-2);

	DM.Set_Mass(mDM_original);
	DM.Set_Interaction_Parameter(coupling_original, detector.Target_Particles());
	return exp(log_limit);
}

void Solar_Reflection_Limit::Compute_Limit_Curve(std::string& ID, obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank)
{
	std::ofstream f;
	if(mpi_rank == 0)
	{
		int CL = std::round(100.0 * certainty_level);
		f.open(TOP_LEVEL_DIR "results/" + ID + "/Reflection_Limit_" + std::to_string(CL) + ".txt");
	}
	for(auto& mass : masses)
	{
		double limit = Upper_Limit(mass, DM, detector, solar_model, halo_model, mpi_rank);
		limits.push_back(limit);
		if(mpi_rank == 0)
		{
			std::cout << mass << "\t" << In_Units(limit, cm * cm) << std::endl;
			f << mass << "\t" << In_Units(limit, cm * cm) << std::endl;
		}
	}
	f.close();
}

void Solar_Reflection_Limit::Export_Curve(std::string& ID, int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::vector<std::vector<double>> data(masses.size(), std::vector<double>(2, 0.0));
		for(unsigned int i = 0; i < masses.size(); i++)
			data[i] = {masses[i], limits[i]};
		int CL = std::round(100.0 * certainty_level);
		libphysica::Export_Table(TOP_LEVEL_DIR "results/" + ID + "/Reflection_Limit_" + std::to_string(CL) + ".txt", data, {GeV, cm * cm});
	}
}

}	// namespace DaMaSCUS_SUN