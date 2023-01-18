#include "Parameter_Scan.hpp"

#include <algorithm>
#include <libconfig.h++>
#include <mpi.h>
#include <set>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "Dark_Photon.hpp"
#include "Data_Generation.hpp"
#include "Reflection_Spectrum.hpp"

namespace DaMaSCUS_SUN
{

using namespace libconfig;
using namespace libphysica::natural_units;

// 1. Configuration class for input file, which extends the obscura::Configuration class.

Configuration::Configuration(std::string cfg_filename, int MPI_rank)
{
	cfg_file	 = cfg_filename;
	results_path = "./";

	// 1. Read the cfg file.
	Read_Config_File();

	// 2. Find the run ID, create a folder and copy the cfg file.
	Initialize_Result_Folder(MPI_rank);

	// 3. DM particle
	Construct_DM_Particle();

	// 4. DM Distribution
	Construct_DM_Distribution();

	// 5. DM-detection experiment
	Construct_DM_Detector();

	// 6. Computation of exclusion limits
	Initialize_Parameters();

	// 7. DaMaSCUS specific parameters
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
		run_mode = config.lookup("run_mode").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'run_mode' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	try
	{
		isoreflection_rings = config.lookup("isoreflection_rings");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'isoreflection_rings' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	try
	{
		interpolation_points = config.lookup("interpolation_points");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'interpolation_points' setting in configuration file." << std::endl;
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

	try
	{
		use_medium_effects		   = config.lookup("use_medium_effects");
		std::string DM_form_factor = config.lookup("DM_form_factor").c_str();
		if(!use_medium_effects && DM_form_factor == "Long-Range")
			use_medium_effects = true;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'use_medium_effects' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		std::string DM_form_factor = config.lookup("DM_form_factor").c_str();
		if(DM_form_factor == "Long-Range")
			zeta = config.lookup("zeta");
		else
			zeta = 0.0;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'zeta' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		perform_full_scan = config.lookup("perform_full_scan");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'perform_full_scan' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(run_mode != "Parameter point" && run_mode != "Parameter scan" && run_mode != "Custom")
	{
		std::cerr << "Error in Configuration::Import_Parameter_Scan_Parameter(): Run mode " << run_mode << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void Configuration::Construct_DM_Particle()
{
	double DM_mass, DM_spin, DM_fraction;
	bool DM_light;
	// 3.1 General properties
	try
	{
		DM_mass = config.lookup("DM_mass");
		DM_mass *= MeV;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_mass' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	try
	{
		DM_spin = config.lookup("DM_spin");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_spin' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	try
	{
		DM_fraction = config.lookup("DM_fraction");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_fraction' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	try
	{
		DM_light = config.lookup("DM_light");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_light' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// 3.2 DM interactions
	std::string DM_interaction;
	try
	{
		DM_interaction = config.lookup("DM_interaction").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_interaction' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// 3.2.1 SI and SD
	if(DM_interaction == "SI" || DM_interaction == "SD")
		Configuration::Construct_DM_Particle_Standard(DM_interaction);
	else if(DM_interaction == "Dark photon")
		Configuration::Construct_DM_Particle_Dark_Photon();
	else
	{
		std::cerr << "Error in DaMaSCUS_SUN::Configuration::Construct_DM_Particle(): 'DM_interaction' setting " << DM_interaction << " in configuration file not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	DM->Set_Mass(DM_mass);
	DM->Set_Spin(DM_spin);
	DM->Set_Fractional_Density(DM_fraction);
	DM->Set_Low_Mass_Mode(DM_light);
}

void Configuration::Construct_DM_Particle_Dark_Photon()
{
	DM = new DM_Particle_Dark_Photon();

	// DM form factor
	std::string DM_form_factor;
	double DM_mediator_mass = -1.0;
	try
	{
		DM_form_factor = config.lookup("DM_form_factor").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "Error in DaMaSCUS_SUN::Configuration::Construct_DM_Particle_DP(): No 'DM_form_factor' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(DM_form_factor == "General")
	{
		try
		{
			DM_mediator_mass = config.lookup("DM_mediator_mass");
			DM_mediator_mass *= MeV;
			if(DM_mediator_mass < 1e-60)
			{
				std::cout << "Error in Configuration::Construct_DM_Particle_Dark_Photon:\tMediator mass for \"General\" form factor needs to be positive. (DM_mediator_mass = " << DM_mediator_mass / MeV << " MeV)\n"
						  << "\tFor long range interactions or zero mediator mass, use \"Long range\" form factor." << std::endl;
				std::exit(EXIT_FAILURE);
			}
		}
		catch(const SettingNotFoundException& nfex)
		{
			std::cerr << "Error in Configuration::Construct_DM_Particle_DP(): No 'DM_mediator_mass' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
	dynamic_cast<DM_Particle_Dark_Photon*>(DM)->Set_FormFactor_DM(DM_form_factor, DM_mediator_mass);

	double DM_cross_section_electron;
	try
	{
		DM_cross_section_electron = config.lookup("DM_cross_section_electron");
		DM_cross_section_electron *= cm * cm;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "Error in Configuration::Construct_DM_Particle_DP(): No 'DM_cross_section_electron' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	DM->Set_Interaction_Parameter(DM_cross_section_electron, "Electrons");
}

void Configuration::Print_Summary(int mpi_rank)
{
	if(mpi_rank == 0)
	{
		Print_Summary_Base(mpi_rank);
		std::cout << "DaMaSCUS-SUN options" << std::endl
				  << std::endl
				  << "\tRun mode:\t\t\t" << run_mode << std::endl
				  << "\tSample size:\t\t\t" << sample_size << std::endl
				  << "\tSc. rate interpolation:\t\t" << ((interpolation_points > 0) ? "[x] (Grid: " + std::to_string(interpolation_points) + "×" + std::to_string(interpolation_points) + ")" : "[ ]") << std::endl
				  << "\tMedium effects:\t\t\t" << (use_medium_effects ? "[x]" : "[ ]") << std::endl;
		if(zeta > 0.0)
			std::cout << "\tQ-cutoff parameter zeta:\t" << zeta << std::endl;
		if(run_mode == "Parameter point" && isoreflection_rings > 1)
			std::cout << "\tIsoreflection rings:\t\t" << isoreflection_rings << std::endl;
		else if(run_mode == "Parameter scan")
			std::cout
				<< "\tCross section (min) [cm^2]:\t" << libphysica::Round(In_Units(cross_section_min, cm * cm)) << std::endl
				<< "\tCross section (max) [cm^2]:\t" << libphysica::Round(In_Units(cross_section_max, cm * cm)) << std::endl
				<< "\tCross section steps:\t\t" << cross_sections << std::endl;
		std::cout << SEPARATOR << std::endl;
	}
}

double Compute_p_Value(unsigned int sample_size, obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, unsigned int rate_interpolation_points, int mpi_rank)
{
	double u_min = detector.Minimum_DM_Speed(DM);

	solar_model.Interpolate_Total_DM_Scattering_Rate(DM, rate_interpolation_points, rate_interpolation_points);
	Simulation_Data data_set(sample_size, u_min);
	data_set.Generate_Data(DM, solar_model, halo_model);
	data_set.Print_Summary(mpi_rank);
	Reflection_Spectrum spectrum(data_set, solar_model, halo_model, DM.mass);
	double p = detector.P_Value(DM, spectrum);
	return (p < 1.0e-100) ? 0.0 : p;
}

// 2. 	Class to perform parameter scans in the (m_DM, sigma)-plane to search for equal-p-value contours.
Parameter_Scan::Parameter_Scan(const std::vector<double>& masses, const std::vector<double>& coupl, std::string ID, unsigned int samplesize, unsigned int interpolation_points, double CL)
: DM_masses(masses), couplings(coupl), sample_size(samplesize), scattering_rate_interpolation_points(interpolation_points), certainty_level(CL)
{
	results_path = TOP_LEVEL_DIR "results/" + ID + "/";
	p_value_grid = std::vector<std::vector<double>>(couplings.size(), std::vector<double>(DM_masses.size(), -1.0));

	// Try to import previous results from an incomplete run
	Import_P_Values();

	std::sort(DM_masses.begin(), DM_masses.end());
	std::sort(couplings.begin(), couplings.end());
}

Parameter_Scan::Parameter_Scan(Configuration& config)
: Parameter_Scan(libphysica::Log_Space(config.constraints_mass_min, config.constraints_mass_max, config.constraints_masses), libphysica::Log_Space(config.cross_section_min, config.cross_section_max, config.cross_sections), config.ID, config.sample_size, config.interpolation_points, config.constraints_certainty)
{
}

void Parameter_Scan::Import_P_Values()
{
	// Import p-values if a corresponding file exists and the grid dimensions fit.
	// CAREFUL: Changes in the grid's mininum/maximum mass/cross section will not be detected at this point.
	std::string filepath = results_path + "P_Values_Grid.txt";
	if(libphysica::File_Exists(filepath))
	{
		std::vector<std::vector<double>> imported_table = libphysica::Import_Table(filepath);
		if(imported_table.size() == p_value_grid.size() && imported_table[0].size() == p_value_grid[0].size())
			p_value_grid = imported_table;
	}
}

bool Parameter_Scan::STA_Point_On_Grid(int row, int column)
{
	return row >= 0 && column >= 0 && row < couplings.size() && column < DM_masses.size();
}

void Parameter_Scan::STA_Go_Forward(int& row, int& column, std::string& STA_direction)
{
	if(STA_direction == "N")
		row++;
	else if(STA_direction == "E")
		column++;
	else if(STA_direction == "S")
		row--;
	else if(STA_direction == "W")
		column--;
}

void Parameter_Scan::STA_Go_Left(int& row, int& column, std::string& STA_direction)
{
	if(STA_direction == "N")
	{
		column--;
		STA_direction = "W";
	}
	else if(STA_direction == "E")
	{
		row++;
		STA_direction = "N";
	}
	else if(STA_direction == "S")
	{
		column++;
		STA_direction = "E";
	}
	else if(STA_direction == "W")
	{
		row--;
		STA_direction = "S";
	}
}

void Parameter_Scan::STA_Go_Right(int& row, int& column, std::string& STA_direction)
{
	if(STA_direction == "N")
	{
		column++;
		STA_direction = "E";
	}
	else if(STA_direction == "E")
	{
		row--;
		STA_direction = "S";
	}
	else if(STA_direction == "S")
	{
		column--;
		STA_direction = "W";
	}
	else if(STA_direction == "W")
	{
		row++;
		STA_direction = "N";
	}
}

void Parameter_Scan::STA_Fill_Gaps()
{
	for(unsigned int row = 0; row < couplings.size(); row++)
	{
		bool excluded = false;
		for(unsigned int column = 0; column < DM_masses.size(); column++)
		{
			double p = p_value_grid[row][column];
			if(p < 0)
				p_value_grid[row][column] = excluded ? 0.0 : 1.0;
			else if(p < 1.0 - certainty_level)
				excluded = true;
			else
				excluded = false;
		}
	}
}

std::vector<double> Parameter_Scan::Find_Contour_Point(int row, int column, int row_previous, int column_previous, double p_critical)
{
	if(row_previous == row)
	{
		double sigma = couplings[row];
		int column_1 = (DM_masses[column] < DM_masses[column_previous]) ? column : column_previous;
		int column_2 = (column_1 == column) ? column_previous : column;
		double x_1	 = DM_masses[column_1];
		double x_2	 = DM_masses[column_2];
		double p_1	 = p_value_grid[row][column_1];
		double p_2	 = p_value_grid[row][column_2];
		double y_1	 = p_1 < 1.0e-100 ? -100.0 : log10(p_1);
		double y_2	 = p_2 < 1.0e-100 ? -100.0 : log10(p_2);
		double y	 = log10(p_critical);
		double x	 = (x_2 - x_1) * (y - y_1) / (y_2 - y_1) + x_1;
		return {x, sigma};
	}
	else
	{
		double mDM = DM_masses[column];
		int row_1  = (couplings[row] < couplings[row_previous]) ? row : row_previous;
		int row_2  = (row_1 == row) ? row_previous : row;
		double x_1 = couplings[row_1];
		double x_2 = couplings[row_2];
		double p_1 = p_value_grid[row_1][column];
		double p_2 = p_value_grid[row_2][column];
		double y_1 = p_1 < 1.0e-100 ? -100.0 : log10(p_1);
		double y_2 = p_2 < 1.0e-100 ? -100.0 : log10(p_2);
		double y   = log10(p_critical);
		double x   = (x_2 - x_1) * (y - y_1) / (y_2 - y_1) + x_1;
		return {mDM, x};
	}
}

std::vector<std::vector<double>> Parameter_Scan::Limit_Curve()
{
	std::vector<std::vector<double>> limit_curve;
	std::string STA_direction = "W";
	int row					  = couplings.size() - 1;
	int column				  = DM_masses.size() - 1;
	int row_previous = -10, column_previous = -10;
	double p_previous = 10.0;
	double p_critical = 1.0 - certainty_level;
	std::vector<int> first_excluded_point;
	unsigned int first_excluded_point_visits = 0;
	while(first_excluded_point_visits < 2)
	{
		double p = STA_Point_On_Grid(row, column) ? p_value_grid[row][column] : 1.0;
		// Abort if no point in the upper row can be excluded:
		if(first_excluded_point.empty() && p > p_critical && row == couplings.size() - 1 && column == 0)
			break;
		// Save the first excluded point and count how often we re-visit that point
		if(p < p_critical && first_excluded_point.empty())
			first_excluded_point = {row, column};
		if(!first_excluded_point.empty() && row == first_excluded_point[0] && column == first_excluded_point[1])
			first_excluded_point_visits++;
		// Interpolate at the boundary to find the point where p == p_critical
		if(STA_Point_On_Grid(row, column) && STA_Point_On_Grid(row_previous, column_previous) && (p - p_critical) * (p_previous - p_critical) < 0.0)
		{
			std::vector<double> contour_point = Find_Contour_Point(row, column, row_previous, column_previous, p_critical);
			if(limit_curve.empty() || limit_curve.back()[1] != couplings[row] || limit_curve.back()[0] != DM_masses[column])
				limit_curve.push_back(contour_point);
		}
		p_previous		= p;
		row_previous	= row;
		column_previous = column;
		if(first_excluded_point.empty())
			STA_Go_Forward(row, column, STA_direction);
		else if(p < p_critical)
			STA_Go_Left(row, column, STA_direction);
		else
			STA_Go_Right(row, column, STA_direction);
	}
	std::reverse(limit_curve.begin(), limit_curve.end());
	return limit_curve;
}

void Parameter_Scan::Perform_STA_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank)
{
	Import_P_Values();
	double mDM_original		 = DM.mass;
	double coupling_original = DM.Get_Interaction_Parameter(detector.Target_Particles());

	std::vector<int> first_excluded_point;

	int counter				  = 0;
	std::string STA_direction = "W";
	int row					  = couplings.size() - 1;
	int column				  = DM_masses.size() - 1;
	double p_critical		  = 1.0 - certainty_level;

	int first_excluded_point_counter = 0;
	while(first_excluded_point_counter < 2)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		double p;
		if(!STA_Point_On_Grid(row, column))
			p = 1.0;
		else if(p_value_grid[row][column] >= 0)
			p = p_value_grid[row][column];
		else
		{
			DM.Set_Interaction_Parameter(couplings[row], detector.Target_Particles());
			DM.Set_Mass(DM_masses[column]);
			double u_min = detector.Minimum_DM_Speed(DM);
			if(mpi_rank == 0)
				std::cout << std::endl
						  << ++counter << ")\t"
						  << "m_DM [MeV]:\t" << libphysica::Round(In_Units(DM.mass, MeV)) << "\t\t"
						  << "sigma_p [cm2]:\t" << libphysica::Round(In_Units(DM.Get_Interaction_Parameter("Nuclei"), cm * cm)) << std::endl
						  << "\tu_min [km/sec]:\t" << libphysica::Round(In_Units(u_min, km / sec)) << "\t\t"
						  << "sigma_e [cm2]:\t" << libphysica::Round(In_Units(DM.Get_Interaction_Parameter("Electrons"), cm * cm)) << std::endl
						  << std::endl;
			Print_Grid(mpi_rank, row, column);
			MPI_Barrier(MPI_COMM_WORLD);

			p = Compute_p_Value(sample_size, DM, detector, solar_model, halo_model, scattering_rate_interpolation_points, mpi_rank);

			p_value_grid[row][column] = p;
			libphysica::Export_Table(results_path + "P_Values_Grid.txt", p_value_grid);
			if(mpi_rank == 0)
			{
				std::cout << std::endl
						  << std::endl;
				libphysica::Print_Box("p = " + std::to_string(libphysica::Round(p)), 1);
			}
		}
		// If the upper row does not contain excluded points, we abort.
		if(first_excluded_point.empty() && p > p_critical && row == couplings.size() - 1 && column == 0)
			break;
		// Check if we arrived back at the first excluded point
		if(first_excluded_point.empty() && p < p_critical)
			first_excluded_point = {row, column};
		if(!first_excluded_point.empty() && row == first_excluded_point[0] && column == first_excluded_point[1])
			first_excluded_point_counter++;

		// Go to the next parameter point
		if(first_excluded_point.empty())
			STA_Go_Forward(row, column, STA_direction);
		else if(p < p_critical)
			STA_Go_Left(row, column, STA_direction);
		else
			STA_Go_Right(row, column, STA_direction);
	}
	STA_Fill_Gaps();
	Print_Grid(mpi_rank);
	libphysica::Export_Table(results_path + "P_Values_Grid.txt", p_value_grid);
	DM.Set_Mass(mDM_original);
	DM.Set_Interaction_Parameter(coupling_original, detector.Target_Particles());
}

void Parameter_Scan::Perform_Full_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank)
{
	Import_P_Values();
	double mDM_original		 = DM.mass;
	double coupling_original = DM.Get_Interaction_Parameter(detector.Target_Particles());

	double p_critical					  = 1.0 - certainty_level;
	unsigned int counter				  = 0;
	unsigned int last_excluded_mass_index = DM_masses.size();
	for(unsigned int i = 0; i < couplings.size(); i++)
	{
		int row			   = couplings.size() - 1 - i;
		bool row_exclusion = false;
		for(unsigned int j = 0; j < DM_masses.size(); j++)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			int column = DM_masses.size() - 1 - j;
			DM.Set_Mass(DM_masses[column]);
			DM.Set_Interaction_Parameter(couplings[row], detector.Target_Particles());
			double u_min = detector.Minimum_DM_Speed(DM);
			double p;
			if(p_value_grid[row][column] >= 0)
				p = p_value_grid[row][column];
			else
			{
				if(mpi_rank == 0)
					std::cout << std::endl
							  << ++counter << ")\t"
							  << "m_DM [MeV]:\t" << libphysica::Round(In_Units(DM.mass, MeV)) << "\t\t"
							  << "sigma_p [cm2]:\t" << libphysica::Round(In_Units(DM.Get_Interaction_Parameter("Nuclei"), cm * cm)) << std::endl
							  << "\tu_min [km/sec]:\t" << libphysica::Round(In_Units(u_min, km / sec)) << "\t\t"
							  << "sigma_e [cm2]:\t" << libphysica::Round(In_Units(DM.Get_Interaction_Parameter("Electrons"), cm * cm)) << std::endl
							  << std::endl;
				Print_Grid(mpi_rank, row, column);
				MPI_Barrier(MPI_COMM_WORLD);

				p = Compute_p_Value(sample_size, DM, detector, solar_model, halo_model, scattering_rate_interpolation_points, mpi_rank);

				p_value_grid[row][column] = p;
				libphysica::Export_Table(results_path + "P_Values_Grid.txt", p_value_grid);
				if(mpi_rank == 0)
				{
					std::cout << std::endl
							  << std::endl;
					libphysica::Print_Box("p = " + std::to_string(libphysica::Round(p)), 1);
				}
			}

			if(p < p_critical)
			{
				row_exclusion			 = true;
				last_excluded_mass_index = j;
			}
		}
		if(!row_exclusion)
		{
			for(int k = 0; k < row; k++)
				for(unsigned int j = 0; j < DM_masses.size(); j++)
					p_value_grid[k][j] = 1.0;
			break;
		}
	}
	libphysica::Export_Table(results_path + "P_Values_Grid.txt", p_value_grid);

	DM.Set_Mass(mDM_original);
	DM.Set_Interaction_Parameter(coupling_original, detector.Target_Particles());
}

void Parameter_Scan::Print_Grid(int mpi_rank, int marker_row, int marker_column)
{
	if(mpi_rank == 0)
	{
		double p_critical = 1.0 - certainty_level;
		std::cout << "\t┌";
		for(unsigned int column = 0; column < DM_masses.size(); column++)
			std::cout << "─";
		std::cout << "┐ σ [cm^2]" << std::endl;
		for(unsigned int row = 0; row < couplings.size(); row++)
		{
			std::cout << "\t┤";
			for(unsigned int column = 0; column < DM_masses.size(); column++)
			{
				double p = p_value_grid[couplings.size() - 1 - row][column];
				if(row == couplings.size() - 1 - marker_row && column == marker_column)
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
		for(unsigned int column = 0; column < DM_masses.size(); column++)
			std::cout << ((column == 0 || column == DM_masses.size() - 1) ? "┬" : "─");
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

void Parameter_Scan::Export_Results(int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::vector<std::vector<double>> table;
		for(unsigned int i = 0; i < DM_masses.size(); i++)
			for(unsigned int j = 0; j < couplings.size(); j++)
				table.push_back({DM_masses[i], couplings[j], p_value_grid[j][i]});
		libphysica::Export_Table(results_path + "P_Values_List.txt", table, {GeV, cm * cm, 1.0});
		int CL										   = std::round(100.0 * certainty_level);
		std::vector<std::vector<double>> limit_contour = Limit_Curve();
		libphysica::Export_Table(results_path + "Reflection_Limit_" + std::to_string(CL) + ".txt", limit_contour, {GeV, cm * cm});
	}
}

}	// namespace DaMaSCUS_SUN