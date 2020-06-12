#include "Data_Generation.hpp"

#include <algorithm>

// Headers from libphysica
#include "Natural_Units.hpp"

using namespace libphysica::natural_units;

Simulation_Data_Set::Simulation_Data_Set()
: number_of_trajectories(0), number_of_free_particles(0), number_of_reflected_particles(0), number_of_captured_particles(0), simulation_time(0.0)
{
}

unsigned int Simulation_Data_Set::Sampe_Size(unsigned int isoreflection_ring)
{
	return data[isoreflection_ring].size();
}

double Simulation_Data_Set::Free_Ratio()
{
	return 1.0 * number_of_free_particles / number_of_trajectories;
}
double Simulation_Data_Set::Capture_Ratio()
{
	return 1.0 * number_of_captured_particles / number_of_trajectories;
}
double Simulation_Data_Set::Reflection_Ratio()
{
	return 1.0 * number_of_reflected_particles / number_of_trajectories;
}

void Simulation_Data_Set::Print_Summary(unsigned int MPI_rank)
{
}

Simulation_Data_Set Generate_Data(unsigned int sample_size, obscura::DM_Particle& DM, Solar_Model& solar_model, double u_min, unsigned int isoreflection_rings)
{
	Trajectory_Simulator simulator(solar_model);
	Simulation_Data_Set data_set;

	obscura::Standard_Halo_Model SHM;
	SHM.Set_Observer_Velocity(libphysica::Vector({0, 0, 0}));

	std::vector<unsigned long int> data_count(isoreflection_rings, 0);
	while(*std::min_element(std::begin(data_count), std::end(data_count)) < sample_size)
	{
		Event IC = Initial_Conditions(SHM, solar_model, simulator.PRNG);
		Hyperbolic_Kepler_Shift(IC, 1.5 * rSun);

		Trajectory_Result result = simulator.Simulate(IC, DM);
	}
	return Simulation_Data_Set();
}
