#include "Reflection_Spectrum.hpp"

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Statistics.hpp"

using namespace libphysica::natural_units;

Reflection_Spectrum::Reflection_Spectrum(const Simulation_Data& simulation_data, unsigned int iso_ring)
: DM_Distribution("Reflection spectrum", 0.0, simulation_data.Lowest_Speed(iso_ring), 1.05 * simulation_data.Highest_Speed(iso_ring))
{
	spectrum_interpolation = libphysica::Perform_KDE(simulation_data.data[iso_ring], v_domain[0], v_domain[1]);
}

double Reflection_Spectrum::PDF_Speed(double v)
{
	return 0.0;
}

double Reflection_Spectrum::Differential_Spectrum(double u)
{
	if(u < v_domain[0] || u > v_domain[1])
		return 0.0;
	else
		return spectrum_interpolation(u);
}

double Reflection_Spectrum::Differential_Flux(double u)
{
	double distance = 1.0 * AU;
	return 1.0 / 4.0 / M_PI / distance / distance * Differential_Spectrum(u);
}

void Reflection_Spectrum::Print_Summary(int mpi_rank)
{
	std::cout << SEPARATOR;
	Print_Summary_Base(mpi_rank);
	std::cout << SEPARATOR;
}

double DM_Entering_Rate(Solar_Model& solar_model, obscura::DM_Distribution& halo_model, double mDM)
{
	double number_density = halo_model.DM_density / mDM;
	double u_average	  = halo_model.Average_Speed();
	double u_inv_average  = halo_model.Eta_Function(0.0);
	double v_esc		  = solar_model.Local_Escape_Speed(rSun);
	return rSun * rSun * M_PI * number_density * (u_average + v_esc * v_esc * u_inv_average);
}
