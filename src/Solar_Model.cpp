#include "Solar_Model.hpp"

#include <cmath>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Statistics.hpp"
#include "Utilities.hpp"

#include "version.hpp"

using namespace libphysica::natural_units;
// 1. Nuclear targets in the Sun
Solar_Isotope::Solar_Isotope(const obscura::Isotope& isotope, const std::vector<std::vector<double>>& density_table, double abundance)
: Isotope(isotope), number_density(libphysica::Interpolation(density_table))
{
	number_density.Multiply(abundance);
}

double Solar_Isotope::Number_Density(double r)
{
	if(r > rSun)
		return 0.0;
	else
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
	std::vector<std::vector<double>> table_vesc(raw_data.size(), std::vector<double>(2, 0.0));
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

std::vector<std::vector<double>> Solar_Model::Create_Number_Density_Table_Electron()
{
	std::vector<std::vector<double>> table(raw_data.size(), std::vector<double>(2, 0.0));
	for(auto& isotope : target_isotopes)
	{
		for(unsigned int i = 0; i < raw_data.size(); i++)
		{
			double r	= raw_data[i][1];
			double n	= isotope.Number_Density(r);
			table[i][0] = r;
			table[i][1] += isotope.Z * n;
		}
	}
	return table;
}

Solar_Model::Solar_Model()
: using_interpolated_rate(false), name("Standard Solar Model AGSS09")
{
	Import_Raw_Data();

	// Interpolate tables.
	mass					   = libphysica::Interpolation(Create_Interpolation_Table(1));
	temperature				   = libphysica::Interpolation(Create_Interpolation_Table(3));
	local_escape_speed_squared = libphysica::Interpolation(Create_Escape_Speed_Table());

	// Nuclear abundances
	obscura::Import_Nuclear_Data();
	std::vector<int> Zs						   = {1, 2, 2, 6, 6, 7, 7, 8, 8, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
	std::vector<double> As					   = {1.0, 4.0, 3.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 23.0, 24.0, 27.0, 28.0, 31.0, 32.0, 35.0, 40.0, 39.0, 40.0, 45.0, 48.0, 51.0, 52.0, 55.0, 56.0, 59.0, 58.0};
	std::vector<unsigned int> included_targets = {0, 1, 7, 26};
	for(auto& target_index : included_targets)
	{
		int Z = Zs[target_index];
		if(target_index < 10)
		{
			double A				 = As[target_index];
			obscura::Isotope isotope = obscura::Get_Isotope(Z, A);
			isotope.abundance		 = 1.0;
			target_isotopes.push_back(Solar_Isotope(isotope, Create_Number_Density_Table(target_index, isotope.mass)));
		}
		else
		{
			obscura::Element element = obscura::Get_Element(Z);
			for(const auto& isotope : element.isotopes)
				target_isotopes.push_back(Solar_Isotope(isotope, Create_Number_Density_Table(target_index, isotope.mass), isotope.abundance));
		}
	}
	// Electron number density
	number_density_electron = libphysica::Interpolation(Create_Number_Density_Table_Electron());
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

double Solar_Model::Number_Density_Nucleus(double r, unsigned int nucleus_index)
{
	if(nucleus_index >= target_isotopes.size())
	{
		std::cerr << "Error in Solar_Model::Number_Density_Nucleus(): Index = " << nucleus_index << " is out of bound (number of targets: " << target_isotopes.size() << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
		return target_isotopes[nucleus_index].Number_Density(r);
}

double Solar_Model::Number_Density_Electron(double r)
{
	if(r > rSun)
		return 0.0;
	else
		return number_density_electron(r);
}

double Solar_Model::DM_Scattering_Rate_Electron(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(r > rSun)
		return 0.0;
	else
	{
		double v_rel = Thermal_Averaged_Relative_Speed(Temperature(r), mElectron, DM.mass);
		return Number_Density_Electron(r) * DM.Sigma_Electron() * v_rel;
	}
}

double Solar_Model::DM_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double r, double DM_speed, unsigned int nucleus_index)
{
	if(nucleus_index >= target_isotopes.size())
	{
		std::cerr << "Error in Solar_Model::Number_Density_Nucleus(): Index = " << nucleus_index << " is out of bound (number of targets: " << target_isotopes.size() << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(r > rSun)
		return 0.0;
	else
	{
		double m_target = target_isotopes[nucleus_index].mass;
		double v_rel	= Thermal_Averaged_Relative_Speed(Temperature(r), m_target, DM.mass);
		double vDM		= 1.0e-3;	//NEEDS TO BE FIXED FOR LIGHT MEDIATORS
		return Number_Density_Nucleus(r, nucleus_index) * DM.Sigma_Nucleus(target_isotopes[nucleus_index], vDM) * v_rel;
	}
}

double Solar_Model::Total_DM_Scattering_Rate(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(using_interpolated_rate && DM_speed < rate_interpolation.domain[1][1])
		return Total_DM_Scattering_Rate_Interpolated(DM, r, DM_speed);
	else
		return Total_DM_Scattering_Rate_Computed(DM, r, DM_speed);
}

double Solar_Model::Total_DM_Scattering_Rate_Computed(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(r > rSun)
		return 0.0;
	else
	{
		double total_rate = DM_Scattering_Rate_Electron(DM, r, DM_speed);
		for(unsigned int i = 0; i < target_isotopes.size(); i++)
			total_rate += DM_Scattering_Rate_Nucleus(DM, r, DM_speed, i);
		return total_rate;
	}
}

double Solar_Model::Total_DM_Scattering_Rate_Interpolated(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(r > rSun)
		return 0.0;
	else
		return rate_interpolation(r, DM_speed);
}

void Solar_Model::Interpolate_Total_DM_Scattering_Rate(obscura::DM_Particle& DM, unsigned int N_radius, unsigned N_speed)
{
	if(N_radius == 0 || N_speed == 0)
		using_interpolated_rate = false;
	else
	{
		using_interpolated_rate = true;

		std::vector<double> radii  = libphysica::Linear_Space(0, rSun, N_radius);
		std::vector<double> speeds = libphysica::Linear_Space(0, 0.3, N_speed);
		std::vector<std::vector<double>> rates;
		for(auto& radius : radii)
			for(auto& speed : speeds)
				rates.push_back({radius, speed, Total_DM_Scattering_Rate_Computed(DM, radius, speed)});

		rate_interpolation = libphysica::Interpolation_2D(rates);
	}
}

void Solar_Model::Print_Summary(int mpi_rank) const
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR
				  << "Solar model:\t\t" << name << std::endl
				  << "Nuclear targets:\t" << target_isotopes.size() << std::endl
				  << std::endl
				  << "Isotope\tZ\tA\tAbund.[%]\tSpin\t<sp>\t<sn>"
				  << SEPARATOR_LINE;
		for(auto& isotope : target_isotopes)
			isotope.Print_Summary(mpi_rank);
		std::cout << SEPARATOR;
	}
}

double Thermal_Averaged_Relative_Speed(double temperature, double mass_target, double v_DM)
{
	double kappa		  = sqrt(mass_target / 2.0 / temperature);
	double relative_speed = (1.0 + 2.0 * pow(kappa * v_DM, 2.0)) * erf(kappa * v_DM) / 2.0 / kappa / kappa / v_DM + exp(-pow(kappa * v_DM, 2.0)) / sqrt(M_PI) / kappa;
	return relative_speed;
}