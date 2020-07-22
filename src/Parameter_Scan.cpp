#include "Parameter_Scan.hpp"

#include <algorithm>
#include <libconfig.h++>
#include <mpi.h>
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

// 1. Configuration class for input file, which extends the obscura::Configuration class.

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

// 2. 	Class to perform parameter scans in the (m_DM, sigma)-plane to search for equal-p-value contours.
Parameter_Scan::Parameter_Scan(Configuration& config)
: Parameter_Scan(libphysica::Log_Space(config.constraints_mass_min, config.constraints_mass_max, config.constraints_masses), libphysica::Log_Space(config.cross_section_min, config.cross_section_max, config.cross_sections), config.sample_size, config.constraints_certainty)
{
}

Parameter_Scan::Parameter_Scan(const std::vector<double>& masses, const std::vector<double>& coupl, unsigned int samplesize, double CL)
: DM_masses(masses), couplings(coupl), sample_size(samplesize), certainty_level(CL)
{
	std::sort(DM_masses.begin(), DM_masses.end());
	std::sort(couplings.begin(), couplings.end());
	p_value_grid = std::vector<std::vector<double>>(couplings.size(), std::vector<double>(DM_masses.size(), -1.0));
}

void Parameter_Scan::STA_Go_Forward(int& row, int& col, std::string& direction)
{
	if(direction == "N")
	{
		if(row == couplings.size() - 1)
			STA_Go_Right(row, col, direction);
		else
			row++;
	}
	else if(direction == "E")
	{
		if(col == DM_masses.size() - 1)
			STA_Go_Right(row, col, direction);
		else
			col++;
	}
	else if(direction == "S")
	{
		if(row == 0)
			STA_Go_Right(row, col, direction);
		else
			row--;
	}
	else if(direction == "W")
	{
		if(col == 0)
			STA_Go_Right(row, col, direction);
		else
			col--;
	}
}

void Parameter_Scan::STA_Go_Left(int& row, int& col, std::string& direction)
{
	if(direction == "N")
	{
		if(col == 0)
			STA_Go_Forward(row, col, direction);
		else
		{
			col--;
			direction = "W";
		}
	}
	else if(direction == "E")
	{
		if(row == couplings.size() - 1)
			STA_Go_Forward(row, col, direction);
		else
		{
			row++;
			direction = "N";
		}
	}
	else if(direction == "S")
	{
		if(col == DM_masses.size() - 1)
			STA_Go_Forward(row, col, direction);
		else
		{
			col++;
			direction = "E";
		}
	}
	else if(direction == "W")
	{
		if(row == 0)
			STA_Go_Forward(row, col, direction);
		else
		{
			row--;
			direction = "S";
		}
	}
}

void Parameter_Scan::STA_Go_Right(int& row, int& col, std::string& direction)
{
	if(direction == "N")
	{
		if(col == DM_masses.size() - 1)
			STA_Go_Forward(row, col, direction);
		else
		{
			col++;
			direction = "E";
		}
	}
	else if(direction == "E")
	{
		if(row == 0)
			STA_Go_Forward(row, col, direction);
		else
		{
			row--;
			direction = "S";
		}
	}
	else if(direction == "S")
	{
		if(col == 0)
			STA_Go_Forward(row, col, direction);
		else
		{
			col--;
			direction = "W";
		}
	}
	else if(direction == "W")
	{
		if(row == couplings.size() - 1)
			STA_Go_Forward(row, col, direction);
		else
		{
			row++;
			direction = "N";
		}
	}
}

void Parameter_Scan::STA_Interpolate_Two_Points(int row, int col, int row_previous, int col_previous, double p_critical)
{
	if(row_previous == row)
	{
		double sigma = couplings[row];
		int col_1	 = (DM_masses[col] < DM_masses[col_previous]) ? col : col_previous;
		int col_2	 = (col_1 == col) ? col_previous : col;
		double x_1	 = DM_masses[col_1];
		double x_2	 = DM_masses[col_2];
		double y_1	 = log10(p_value_grid[row][col_1]);
		double y_2	 = log10(p_value_grid[row][col_2]);
		double y	 = log10(p_critical);
		double x	 = (x_2 - x_1) * (y - y_1) / (y_2 - y_1) + x_1;
		if(limit_curve.empty() || limit_curve.back()[1] != sigma)
			limit_curve.push_back({x, sigma});
	}
	else
	{
		double mDM = DM_masses[col];
		int row_1  = (couplings[row] < couplings[row_previous]) ? row : row_previous;
		int row_2  = (row_1 == row) ? row_previous : row;
		double x_1 = couplings[row_1];
		double x_2 = couplings[row_2];
		double y_1 = log10(p_value_grid[row_1][col]);
		double y_2 = log10(p_value_grid[row_2][col]);
		double y   = log10(p_critical);
		double x   = (x_2 - x_1) * (y - y_1) / (y_2 - y_1) + x_1;
		if(limit_curve.empty() || limit_curve.back()[0] != mDM)
			limit_curve.push_back({mDM, x});
	}
}

void Parameter_Scan::STA_Fill_Gaps()
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
					  << "sigma_p [cm2]:\t" << libphysica::Round(In_Units(DM.Get_Interaction_Parameter("Nuclei"), cm * cm)) << std::endl
					  << "\tu_min [km/sec]:\t" << libphysica::Round(In_Units(u_min, km / sec)) << "\t\t"
					  << "sigma_e [cm2]:\t" << libphysica::Round(In_Units(DM.Get_Interaction_Parameter("Electrons"), cm * cm)) << std::endl
					  << std::endl;
		Print_Grid(mpi_rank, row, col);
		MPI_Barrier(MPI_COMM_WORLD);

		solar_model.Interpolate_Total_DM_Scattering_Rate(DM, 1000, 50);
		Simulation_Data data_set(sample_size, u_min);
		data_set.Generate_Data(DM, solar_model, halo_model);
		Reflection_Spectrum spectrum(data_set, solar_model, halo_model, DM.mass);
		double p = detector.P_Value(DM, spectrum);
		p		 = (p < 1.0e-100) ? 0.0 : p;

		DM.Set_Mass(mDM_original);
		DM.Set_Interaction_Parameter(coupling_original, detector.Target_Particles());

		p_value_grid[row][col] = p;
		if(mpi_rank == 0)
		{
			std::cout << std::endl
					  << std::endl;
			libphysica::Print_Box("p = " + std::to_string(libphysica::Round(p)), 1);
		}
		return p;
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
			MPI_Barrier(MPI_COMM_WORLD);
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

void Parameter_Scan::Perform_STA_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank)
{
	std::string direction = "S";
	int row				  = couplings.size() - 1;
	int col				  = DM_masses.size() - 1;
	int row_previous, col_previous;
	double p_previous					  = 10.0;
	double p_critical					  = 1.0 - certainty_level;
	std::vector<int> first_excluded_point = {-1, -1};
	int first_excluded_point_counter	  = 0;

	while(first_excluded_point_counter < 2)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		double p = Compute_p_Value(row, col, DM, detector, solar_model, halo_model, mpi_rank);

		// Stopping criterion
		if(p < p_critical && first_excluded_point[0] < 0)
			first_excluded_point = {row, col};
		if(row == first_excluded_point[0] && col == first_excluded_point[1])
			first_excluded_point_counter++;

		// Interpolate at the boundary to find the point where p == p_critical
		if(counter > 1 && (p - p_critical) * (p_previous - p_critical) < 0.0)
			STA_Interpolate_Two_Points(row, col, row_previous, col_previous, p_critical);
		p_previous	 = p;
		row_previous = row;
		col_previous = col;
		if(p < p_critical)
			STA_Go_Left(row, col, direction);
		else if(col == 0 && direction == "S")
		{
			direction = "N";
			STA_Go_Forward(row, col, direction);
		}
		else if(row == 0 && direction == "E")
		{
			direction = "W";
			STA_Go_Forward(row, col, direction);
		}
		else
			STA_Go_Right(row, col, direction);
	}
	STA_Fill_Gaps();
	std::reverse(limit_curve.begin(), limit_curve.end());
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

void Parameter_Scan::Export_Results(const std::string& ID, int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::vector<std::vector<double>> table;
		for(unsigned int i = 0; i < DM_masses.size(); i++)
			for(unsigned int j = 0; j < couplings.size(); j++)
				table.push_back({DM_masses[i], couplings[j], p_value_grid[j][i]});
		libphysica::Export_Table(TOP_LEVEL_DIR "results/" + ID + "/P_Values_List.txt", table, {GeV, cm * cm, 1.0});
		libphysica::Export_Table(TOP_LEVEL_DIR "results/" + ID + "/P_Values_Grid.txt", p_value_grid);
		int CL = std::round(100.0 * certainty_level);
		libphysica::Export_Table(TOP_LEVEL_DIR "results/" + ID + "/Reflection_Limit_" + std::to_string(CL) + ".txt", limit_curve, {GeV, cm * cm});
	}
}

void Parameter_Scan::Print_Grid(int mpi_rank, int marker_row, int marker_col)
{
	if(mpi_rank == 0)
	{
		double p_critical = 1.0 - certainty_level;
		std::cout << "\t┌";
		for(int col = 0; col < DM_masses.size(); col++)
			std::cout << "─";
		std::cout << "┐ σ [cm^2]" << std::endl;
		for(unsigned int row = 0; row < couplings.size(); row++)
		{
			std::cout << "\t┤";
			for(unsigned int col = 0; col < DM_masses.size(); col++)
			{
				double p = p_value_grid[couplings.size() - 1 - row][col];
				if(row == couplings.size() - 1 - marker_row && col == marker_col)
					std::cout << "¤";
				else if(p < 0.0)
					std::cout << "·";
				else if(p < p_critical)
					std::cout << "█";
				else if(p > p_critical)
					std::cout << "░";
			}
			std::cout << "├";	// ((row == 0 || row == couplings.size() - 1) ? "├" : "│");
			if(row == 0)
				std::cout << " " << libphysica::Round(In_Units(couplings.back(), cm * cm));
			else if(row == couplings.size() - 1)
				std::cout << " " << libphysica::Round(In_Units(couplings.front(), cm * cm));
			std::cout << std::endl;
		}
		std::cout << "\t└";
		for(unsigned int col = 0; col < DM_masses.size(); col++)
			std::cout << ((col == 0 || col == DM_masses.size() - 1) ? "┬" : "─");
		if(DM_masses.size() > 5)
		{
			std::cout << "┘ mDM[MeV]\n\t ";
			for(unsigned int i = 0; i < DM_masses.size() - 1; i++)
				std::cout
					<< " ";
			std::cout << libphysica::Round(In_Units(DM_masses.back(), MeV))
					  << "\r\t " << libphysica::Round(In_Units(DM_masses.front(), MeV)) << std::endl;
		}
		std::cout << std::endl;
	}
}

}	// namespace DaMaSCUS_SUN