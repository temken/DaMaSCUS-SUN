#include "Reflection_Spectrum.hpp"

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Statistics.hpp"

using namespace libphysica::natural_units;

Reflection_Spectrum::Reflection_Spectrum(const Simulation_Data& simulation_data, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, double mDM, int iso_ring)
: DM_Distribution("Reflection spectrum", 0.0, simulation_data.Minimum_Speed(), 1.05 * simulation_data.Highest_Speed(iso_ring)), distance(AU)
{
	kde_speed								   = libphysica::Perform_KDE(simulation_data.data[iso_ring], v_domain[0], v_domain[1]);
	total_entering_rate						   = DM_Entering_Rate(solar_model, halo_model, mDM);
	total_reflection_rate					   = simulation_data.Reflection_Ratio(iso_ring) * total_entering_rate;
	unsigned int number_of_isoreflection_rings = simulation_data.data.size();
	total_reflection_flux					   = total_reflection_rate * 1.0 / 4.0 / M_PI / distance / distance * number_of_isoreflection_rings;
}

double Reflection_Spectrum::PDF_Speed(double u)
{
	if(u < v_domain[0] || u > v_domain[1])
		return 0.0;
	else
		return kde_speed(u);
}

double Reflection_Spectrum::Differential_Spectrum(double u)
{
	if(u < v_domain[0] || u > v_domain[1])
		return 0.0;
	else
		return total_reflection_rate * kde_speed(u);
}

double Reflection_Spectrum::Differential_DM_Flux(double u, double mDM)
{
	if(u < v_domain[0] || u > v_domain[1])
		return 0.0;
	else
		return total_reflection_flux * kde_speed(u);
}

void Reflection_Spectrum::Print_Summary(int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR;
		Print_Summary_Base(mpi_rank);
		std::cout << "Total DM entering rate [1/s]:\t\t" << libphysica::Round(In_Units(total_entering_rate, 1.0 / sec)) << std::endl
				  << "Total reflection rate [1/s]:\t\t" << libphysica::Round(In_Units(total_reflection_rate, 1.0 / sec)) << std::endl
				  << "Distance from Sun [AU]:\t\t\t" << libphysica::Round(In_Units(distance, AU)) << std::endl
				  << "Total reflection flux [1/s/cm^2]:\t" << libphysica::Round(In_Units(total_reflection_flux, 1.0 / sec / cm / cm)) << std::endl
				  << SEPARATOR;
	}
}

void Reflection_Spectrum::Set_Distance(double d)
{
	total_reflection_flux *= distance * distance / d / d;
	distance = d;
}

double DM_Entering_Rate(Solar_Model& solar_model, obscura::DM_Distribution& halo_model, double mDM)
{
	double number_density = halo_model.DM_density / mDM;
	double u_average	  = halo_model.Average_Speed();
	double u_inv_average  = halo_model.Eta_Function(0.0);
	double v_esc		  = solar_model.Local_Escape_Speed(rSun);
	return rSun * rSun * M_PI * number_density * (u_average + v_esc * v_esc * u_inv_average);
}
