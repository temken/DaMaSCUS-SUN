#include "Solar_Model.hpp"

#include <cmath>
#include <mpi.h>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

#include "Scattering_Rates.hpp"
#include "version.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;
using namespace std::complex_literals;

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

// 2. Plasma struct to act as target for DM scatterings and describe in-medium effects
Plasma::Plasma(double T, double ne, std::vector<double>& nn, std::vector<obscura::Isotope>& iso)
: temperature(T), number_density_electrons(ne), number_densities_nuclei(nn), nuclei(iso)
{
	libphysica::Check_For_Error(nuclei.size() != number_densities_nuclei.size(), "Plasma::Plasma", "The number of nuclei and the number of densities must be the same.");
}

std::complex<double> Plasma::Polarization_Tensor_L(double q0, double q)
{
	std::complex<double> PI_L = 0.0;
	// 1. Electron contributions
	PI_L += Polarization_Tensor_Longitudinal(q0, q, temperature, number_density_electrons, mElectron, 1.0);

	// 2. Nuclear contributions
	for(unsigned int i = 0; i < nuclei.size(); i++)
		PI_L += Polarization_Tensor_Longitudinal(q0, q, temperature, number_densities_nuclei[i], nuclei[i].mass, nuclei[i].Z);
	return PI_L;
}

double Plasma::Form_Factor_Medium_Effects(double q0, double q)
{
	std::complex<double> denominator = q * q + Polarization_Tensor_L(q0, q);
	return q * q * q * q / std::norm(denominator);
}

std::complex<double> Plasma_Dispersion_Function(double x)
{
	double expmx2_erfi = 2.0 / std::sqrt(M_PI) * libphysica::Dawson_Integral(x);   // remove exp(-x^2) exp(x^2) to avoid inf value
	return std::sqrt(M_PI) * (1i * std::exp(-x * x) - expmx2_erfi);
}

std::complex<double> Polarization_Tensor_Longitudinal(double q0, double q, double temperature, double number_density, double mass, double Z)
{
	bool use_vlaslov_approximation = false;

	double sigma_MB				= std::sqrt(temperature / mass);
	double plasma_frequency_sqr = Elementary_Charge * Elementary_Charge * Z * Z * number_density / mass;
	double xi					= q0 / std::sqrt(2.0) / sigma_MB / q;
	double delta				= q / 2.0 / std::sqrt(2.0) / mass / sigma_MB;

	if(!use_vlaslov_approximation)
		return plasma_frequency_sqr * mass / q0 * xi * (Plasma_Dispersion_Function(xi - delta) - Plasma_Dispersion_Function(xi + delta));
	else
		return plasma_frequency_sqr / sigma_MB / sigma_MB * (1.0 + xi * Plasma_Dispersion_Function(xi));
}

// 3. Solar model
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
		if(i == 0)
			last_line[i] = mSun;
		else if(i == 1)
			last_line[i] = rSun;
		else if(i == 2)
			last_line[i] = 5800 * Kelvin;	// temperature of the photosphere (http://solar-center.stanford.edu/vitalstats.html)
		else if(i == 3)
			last_line[i] = 1.0e-9 * gram / cm / cm / cm;   // mass density of the photosphere (http://solar-center.stanford.edu/vitalstats.html)
		else
			last_line[i] = raw_data.back()[i];
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
		double integral	 = (r < rSun) ? libphysica::Integrate(integrand, r, rSun) : mSun * (1.0 / r - 1.0 / rSun);
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
	for(auto& isotope : all_isotopes)
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

Solar_Model::Solar_Model(bool medium_effects, double q_cutoff_parameter)
: using_interpolated_rate(false), use_medium_effects(medium_effects), zeta(q_cutoff_parameter), name("Standard Solar Model AGSS09")
{
	Import_Raw_Data();

	// Interpolate tables.
	mass					   = libphysica::Interpolation(Create_Interpolation_Table(1));
	temperature				   = libphysica::Interpolation(Create_Interpolation_Table(3));
	local_escape_speed_squared = libphysica::Interpolation(Create_Escape_Speed_Table());
	mass_density			   = libphysica::Interpolation(Create_Interpolation_Table(4));

	// Nuclear abundances
	obscura::Import_Nuclear_Data();
	std::vector<int> Zs	   = {1, 2, 2, 6, 6, 7, 7, 8, 8, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
	std::vector<double> As = {1.0, 4.0, 3.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0};
	for(auto& target_index : libphysica::Range(Zs.size()))
	{
		int Z = Zs[target_index];
		if(target_index < 10)
		{
			double A				 = As[target_index];
			obscura::Isotope isotope = obscura::Get_Isotope(Z, A);
			isotope.abundance		 = 1.0;
			all_isotopes.push_back(Solar_Isotope(isotope, Create_Number_Density_Table(target_index, isotope.mass)));
		}
		else
		{
			obscura::Nucleus nucleus = obscura::Get_Nucleus(Z);
			for(const auto& isotope : nucleus.isotopes)
				all_isotopes.push_back(Solar_Isotope(isotope, Create_Number_Density_Table(target_index, isotope.mass), isotope.abundance));
		}
	}
	// Electron number density
	number_density_electron = libphysica::Interpolation(Create_Number_Density_Table_Electron());
	// Include only a subset of targets
	// std::vector<int> included_target_indices = {};
	std::vector<int> included_target_indices = {0, 1, 2, 7, 54};
	// std::vector<int> included_target_indices = libphysica::Range(all_isotopes.size());
	for(auto& index : included_target_indices)
		target_isotopes.push_back(all_isotopes[index]);
}

double Solar_Model::Mass(double r)
{
	if(r > rSun)
		return mSun;
	else
		return mass(r);
}

double Solar_Model::Mass_Density(double r)
{
	if(r > rSun)
		return 0.0;
	else
		return mass_density(r);
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

double Solar_Model::DM_Scattering_Rate_Electron(obscura::DM_Particle& DM, double r, double vDM)
{
	if(r > rSun)
		return 0.0;
	else
	{
		double electron_density = Number_Density_Electron(r);
		double T				= Temperature(r);
		double qMin				= zeta * DM.mass * vDM;
		double qMax				= 2.0 * libphysica::Reduced_Mass(DM.mass, mElectron) * vRel_max;
		return Total_Scattering_Rate_Electron(DM, vDM, electron_density, T, use_medium_effects, qMin, qMax);
	}
}

double Solar_Model::DM_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double r, double vDM, unsigned int nucleus_index)
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
		double nucleus_density = Number_Density_Nucleus(r, nucleus_index);
		double T			   = Temperature(r);
		double qMin			   = zeta * DM.mass * vDM;
		double qMax			   = 2.0 * libphysica::Reduced_Mass(DM.mass, target_isotopes[nucleus_index].mass) * vRel_max;
		return Total_Scattering_Rate_Nucleus(DM, vDM, target_isotopes[nucleus_index], nucleus_density, T, use_medium_effects, qMin, qMax);
	}
}

double Solar_Model::Total_DM_Scattering_Rate(obscura::DM_Particle& DM, double r, double vDM)
{
	if(using_interpolated_rate && vDM < rate_interpolation.domain[1][1])
		return Total_DM_Scattering_Rate_Interpolated(DM, r, vDM);
	else
	{
		if(using_interpolated_rate)
			std::cerr << "Warning Solar_Model::Total_DM_Scattering_Rate(): DM speed is out of bound (vDM = " << vDM << ")\n\tScattering rate must be computed on the fly." << std::endl;
		return Total_DM_Scattering_Rate_Computed(DM, r, vDM);
	}
}

double Solar_Model::Total_DM_Scattering_Rate_Computed(obscura::DM_Particle& DM, double r, double vDM)
{
	if(r > rSun)
		return 0.0;
	else
	{
		double total_rate = DM_Scattering_Rate_Electron(DM, r, vDM);
		for(unsigned int i = 0; i < target_isotopes.size(); i++)
			total_rate += DM_Scattering_Rate_Nucleus(DM, r, vDM, i);
		return total_rate;
	}
}

double Solar_Model::Total_DM_Scattering_Rate_Interpolated(obscura::DM_Particle& DM, double r, double vDM)
{
	if(r > rSun)
		return 0.0;
	else
		return rate_interpolation(r, vDM);
}

void Solar_Model::Interpolate_Total_DM_Scattering_Rate(obscura::DM_Particle& DM, unsigned int N_radius, unsigned int N_speed)
{
	if(N_radius == 0 || N_speed == 0)
		using_interpolated_rate = false;
	else
	{
		int mpi_processes, mpi_rank;
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

		using_interpolated_rate = true;

		double vMax					 = 0.75;
		unsigned int local_N_radius	 = std::ceil(1.0 * N_radius / mpi_processes);
		unsigned int global_N_radius = mpi_processes * local_N_radius;

		std::vector<double> global_radii = libphysica::Linear_Space(0, rSun, global_N_radius);
		std::vector<double> local_radii(local_N_radius, 0.0);

		// Compute the table in parallel
		MPI_Scatter(global_radii.data(), local_N_radius, MPI_DOUBLE, local_radii.data(), local_N_radius, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		std::vector<double> speeds = libphysica::Linear_Space(0, vMax, N_speed);
		std::vector<double> local_rates;
		std::vector<double> global_rates(N_speed * global_N_radius, 0.0);
		for(auto& radius : local_radii)
			for(auto& speed : speeds)
				local_rates.push_back(Total_DM_Scattering_Rate_Computed(DM, radius, speed));
		MPI_Allgather(local_rates.data(), local_N_radius * N_speed, MPI_DOUBLE, global_rates.data(), local_N_radius * N_speed, MPI_DOUBLE, MPI_COMM_WORLD);

		// Re-organize into a 2D array and interpolate.
		std::vector<std::vector<double>> rates;
		int i = 0;
		for(auto& radius : global_radii)
			for(auto& speed : speeds)
				rates.push_back({radius, speed, global_rates[i++]});
		rate_interpolation = libphysica::Interpolation_2D(rates);
	}
}

void Solar_Model::Print_Summary(int mpi_rank) const
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR
				  << "Solar model:\t\t" << name << std::endl
				  << "In medium effects:\t" << (use_medium_effects ? "[x]" : "[ ]") << std::endl
				  << "Nuclear targets:\t" << target_isotopes.size() << std::endl
				  << std::endl
				  << "Isotope\tZ\tA\tAbund.[%]\tSpin\t<sp>\t<sn>"
				  << SEPARATOR_LINE;
		for(auto& isotope : target_isotopes)
			isotope.Print_Summary(mpi_rank);
		std::cout << SEPARATOR;
	}
}
}	// namespace DaMaSCUS_SUN