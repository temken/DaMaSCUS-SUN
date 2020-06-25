#include "Data_Generation.hpp"

#include <algorithm>
#include <chrono>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Utilities.hpp"

// Headers from obscura
#include "Astronomy.hpp"

using namespace libphysica::natural_units;

Simulation_Data::Simulation_Data(unsigned int sample_size, double u_min, unsigned int iso_rings)
: min_sample_size(sample_size), minimum_speed(u_min), isoreflection_rings(iso_rings), number_of_trajectories(0), number_of_free_particles(0), number_of_reflected_particles(0), number_of_captured_particles(0), average_number_of_scatterings(0.0), computing_time(0.0), number_of_data_points(std::vector<unsigned long int>(iso_rings, 0)), data(iso_rings, std::vector<libphysica::DataPoint>())
{
}

void Simulation_Data::Configure(double initial_radius, unsigned int min_scattering, unsigned int max_scattering, unsigned long int max_free_steps)
{
	initial_and_final_radius	  = initial_radius;
	minimum_number_of_scatterings = min_scattering;
	maximum_number_of_scatterings = max_scattering;
	maximum_free_time_steps		  = max_free_steps;
}

void Simulation_Data::Generate_Data(obscura::DM_Particle& DM, Solar_Model& solar_model, unsigned int fixed_seed)
{
	auto time_start = std::chrono::system_clock::now();

	//Configure the simulator
	Trajectory_Simulator simulator(solar_model, maximum_free_time_steps, maximum_number_of_scatterings, initial_and_final_radius);
	if(fixed_seed != 0)
		simulator.Fix_PRNG_Seed(fixed_seed);

	obscura::Standard_Halo_Model SHM;
	SHM.Set_Observer_Velocity(libphysica::Vector({0, 0, 0}));

	unsigned int smallest_sample_size = 0;
	while(smallest_sample_size < min_sample_size)
	{
		Event IC = Initial_Conditions(SHM, solar_model, simulator.PRNG);
		Hyperbolic_Kepler_Shift(IC, initial_and_final_radius);
		Trajectory_Result trajectory = simulator.Simulate(IC, DM);

		number_of_trajectories++;
		average_number_of_scatterings = 1.0 / number_of_trajectories * ((number_of_trajectories - 1) * average_number_of_scatterings + trajectory.number_of_scatterings);

		if(trajectory.Particle_Captured())
			number_of_captured_particles++;
		else
		{
			if(trajectory.Particle_Free())
				number_of_free_particles++;
			else if(trajectory.Particle_Reflected())
				number_of_reflected_particles++;
			Hyperbolic_Kepler_Shift(trajectory.final_event, 1.0 * AU);
			double v_final = trajectory.final_event.Speed();
			if(trajectory.number_of_scatterings >= minimum_number_of_scatterings && v_final > minimum_speed)
			{
				unsigned int isoreflection_ring = (isoreflection_rings == 1) ? 0 : trajectory.final_event.Isoreflection_Ring(obscura::Sun_Velocity(), isoreflection_rings);
				number_of_data_points[isoreflection_ring]++;
				data[isoreflection_ring].push_back(libphysica::DataPoint(v_final));
				smallest_sample_size = *std::min_element(std::begin(number_of_data_points), std::end(number_of_data_points));
				libphysica::Print_Progress_Bar(1.0 * smallest_sample_size / min_sample_size);
			}
		}
	}
	auto time_end  = std::chrono::system_clock::now();
	computing_time = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
}

double Simulation_Data::Free_Ratio()
{
	return 1.0 * number_of_free_particles / number_of_trajectories;
}
double Simulation_Data::Capture_Ratio()
{
	return 1.0 * number_of_captured_particles / number_of_trajectories;
}
double Simulation_Data::Reflection_Ratio()
{
	return 1.0 * number_of_reflected_particles / number_of_trajectories;
}

void Simulation_Data::Print_Summary(unsigned int MPI_rank)
{
	if(MPI_rank == 0)
	{
		unsigned long int number_of_data_points_tot = std::accumulate(number_of_data_points.begin(), number_of_data_points.end(), 0);
		std::cout << SEPARATOR
				  << "Simulation data summary" << std::endl
				  << std::endl
				  << "Configuration:" << std::endl
				  << "Minimum DM speed [km/sec]:\t" << libphysica::Round(In_Units(minimum_speed, km / sec)) << std::endl
				  << "Minimum sample size:\t\t" << min_sample_size << std::endl
				  << "Isoreflection rings:\t\t" << isoreflection_rings << std::endl
				  << std::endl
				  << "Results:" << std::endl
				  << "Simulated trajectories:\t\t" << number_of_trajectories << std::endl
				  << "Average # of scatterings:\t" << libphysica::Round(average_number_of_scatterings) << std::endl
				  << "Free particles [%]:\t\t" << libphysica::Round(100.0 * Free_Ratio()) << std::endl
				  << "Reflected particles [%]:\t" << libphysica::Round(100.0 * Reflection_Ratio()) << std::endl
				  << "Captured particles [%]:\t\t" << libphysica::Round(100.0 * Capture_Ratio()) << std::endl;

		if(isoreflection_rings > 1)
		{
			std::cout << std::endl
					  << "\tRing\tData points\tRelative [%]\t<u> [km/sec]\tu_max [km/sec]" << std::endl;
			for(unsigned int i = 0; i < isoreflection_rings; i++)
			{
				double rel_number_of_data_points = 100.0 * number_of_data_points[i] / number_of_data_points_tot;
				std::vector<double> u_average	 = libphysica::Weighted_Average(data[i]);
				double u_max					 = (*std::max_element(data[i].begin(), data[i].end())).value;
				std::cout << "\t" << i + 1 << "\t" << number_of_data_points[i] << "\t\t" << libphysica::Round(rel_number_of_data_points) << "\t\t" << libphysica::Round(In_Units(u_average[0], km / sec)) << " +- " << libphysica::Round(In_Units(u_average[1], km / sec)) << "\t" << libphysica::Round(In_Units(u_max, km / sec)) << std::endl;
			}
		}
		else
		{
			std::vector<double> u_average = libphysica::Weighted_Average(data[0]);
			double u_max				  = (*std::max_element(data[0].begin(), data[0].end())).value;
			std::cout << "<u> [km/sec]:\t\t\t" << libphysica::Round(In_Units(u_average[0], km / sec)) << " +- " << libphysica::Round(In_Units(u_average[1], km / sec)) << std::endl
					  << "u_max [km/sec]:\t\t\t" << libphysica::Round(In_Units(u_max, km / sec)) << std::endl;
		}
		std::cout << std::endl;
		if(computing_time > 60.0)

			std::cout << "Simulation time:\t\t"
					  << "[" << floor(computing_time / 3600.0) << ":" << floor(fmod(computing_time / 60.0, 60.0)) << ":" << floor(fmod(computing_time, 60.0)) << "]." << std::endl;
		else

			std::cout << "Simulation time [s]:\t\t" << libphysica::Round(computing_time) << std::endl;
		std::cout << "Trajectory rate [1/s]:\t\t" << libphysica::Round(1.0 * number_of_trajectories / computing_time) << std::endl;
		std::cout << SEPARATOR << std::endl;
	}
}