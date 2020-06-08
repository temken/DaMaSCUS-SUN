#include "Solar_Model.hpp"

#include <cmath>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Statistics.hpp"
#include "Utilities.hpp"

using namespace libphysica::natural_units;
// 1. Nuclear targets in the Sun
Solar_Nucleus::Solar_Nucleus(const obscura::Isotope& isotope, const std::vector<std::vector<double>>& density_table, std::string label)
: Element(isotope), number_density(libphysica::Interpolation(density_table))
{
	if(label != "")
		name = label;
}

Solar_Nucleus::Solar_Nucleus(const std::vector<obscura::Isotope>& isotopes, const std::vector<std::vector<double>>& density_table, std::string label)
: Element(isotopes), number_density(libphysica::Interpolation(density_table))
{
	if(label != "")
		name = label;
}

double Solar_Nucleus::Number_Density(double r)
{
	return number_density(r);
}

// 2. Solar model
// Auxiliary functions for the data import
void Solar_Model::Import_Raw_Data()
{
	std::vector<double> units = {mSun, rSun, Kelvin, gram / cm / cm / cm, dyne / cm / cm, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	raw_data				  = libphysica::Import_Table(PROJECT_DIR "/data/model_agss09.dat", units, 20);
	// Add a first and last line to ensure a full domain of the interpolations.
	std::vector<double> first_line(35, 0.0);
	std::vector<double> last_line(35, 0.0);
	for(unsigned int i = 0; i < 35; i++)
	{
		first_line[i] = (i == 0 || i == 1) ? 0.0 : raw_data.front()[i];
		last_line[i]  = (i == 0 || i == 1) ? ((i == 0) ? mSun : rSun) : raw_data.back()[i];
	}
	raw_data.insert(raw_data.begin(), first_line);
	raw_data.push_back(last_line);
}

std::vector<std::vector<double>> Solar_Model::Create_Interpolation_Table(unsigned int row) const
{
	std::vector<std::vector<double>> table(raw_data.size(), std::vector<double>(2, 0.0));
	for(unsigned int i = 0; i < table.size(); i++)
	{
		table[i][0] = raw_data[i][1];
		table[i][1] = raw_data[i][row - 1];
	}
	return table;
}

std::vector<std::vector<double>> Solar_Model::Create_Escape_Speed_Table()
{
	std::vector<std::vector<double>> table_vesc(raw_data.size(), std::vector<double>(raw_data[0].size(), 0.0));
	for(unsigned int i = 0; i < raw_data.size(); i++)
	{
		double r	   = raw_data[i][1];
		auto integrand = [this](double x) {
			if(x == 0)
				return 0.0;
			else
				return Mass(x) / x / x;
		};
		double integral	 = libphysica::Integrate(integrand, r, rSun, 1.e30);
		table_vesc[i][0] = r;
		table_vesc[i][1] = 2.0 * G_Newton * mSun / rSun * (1.0 + rSun / mSun * integral);
	}
	return table_vesc;
}

std::vector<std::vector<double>> Solar_Model::Create_Number_Density_Table(unsigned int target, double mass) const
{
	std::vector<std::vector<double>> table;
	for(unsigned int i = 0; i < raw_data.size(); i++)
	{
		double r   = raw_data[i][1];
		double rho = raw_data[i][3];
		double f   = raw_data[i][6 + target];
		double n   = f * rho / mass;
		table.push_back({r, n});
	}
	return table;
}

Solar_Model::Solar_Model()
: name("Standard Solar Model AGSS09")
{
	Import_Raw_Data();

	// Interpolate tables.
	mass					   = libphysica::Interpolation(Create_Interpolation_Table(1));
	temperature				   = libphysica::Interpolation(Create_Interpolation_Table(3));
	local_escape_speed_squared = libphysica::Interpolation(Create_Escape_Speed_Table());

	// Nuclear abundances
	obscura::Import_Nuclear_Data();
	std::vector<std::string> names			   = {"H1", "He4", "He3", "C12", "C13", "N14", "N15", "O16", "O17", "O18", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"};
	std::vector<int> Zs						   = {1, 2, 2, 6, 6, 7, 7, 8, 8, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
	std::vector<double> As					   = {1.0, 4.0, 3.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 23.0, 24.0, 27.0, 28.0, 31.0, 32.0, 35.0, 40.0, 39.0, 40.0, 45.0, 48.0, 51.0, 52.0, 55.0, 56.0, 59.0, 58.0};
	std::vector<unsigned int> included_targets = {0, 1, 7, 26};
	for(auto& target : included_targets)
	{
		std::string name = names[target];
		int Z			 = Zs[target];
		if(target < 10)
		{
			double A				 = As[target];
			obscura::Isotope isotope = obscura::Get_Isotope(Z, A);
			isotope.abundance		 = 1.0;
			nuclear_targets.push_back(Solar_Nucleus(isotope, Create_Number_Density_Table(target, isotope.mass), name));
		}
		else
		{
			obscura::Element element = obscura::Get_Element(Z);
			nuclear_targets.push_back(Solar_Nucleus(element.isotopes, Create_Number_Density_Table(target, element.Average_Nuclear_Mass())));
		}
	}
}

double Solar_Model::Mass(double r)
{
	if(r > rSun)
		return mSun;
	else
		return mass(r);
}

double Solar_Model::Temperature(double r)
{
	return temperature(r);
}

double Solar_Model::Local_Escape_Speed(double r)
{
	if(r > rSun)
		return sqrt(2 * G_Newton * mSun / r);
	else
		return sqrt(local_escape_speed_squared(r));
}

void Solar_Model::Print_Summary(int MPI_rank) const
{
	if(MPI_rank == 0)
	{
		std::cout << SEPARATOR
				  << "Solar model:\t\t" << name << std::endl
				  << "Nuclear targets:\t" << nuclear_targets.size();
		// for(auto& target : nuclear_targets)
		// 	target.Print_Summary(MPI_rank);
		std::cout << SEPARATOR;
	}
}