#ifndef __Data_Generation_hpp_
#define __Data_Generation_hpp_

#include <vector>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Statistics.hpp"

// Headers from obscura
#include "DM_Particle.hpp"

#include "Simulation_Trajectory.hpp"

class Simulation_Data
{
  private:
	// Configuration
	unsigned int min_sample_size;
	double minimum_speed;
	unsigned int isoreflection_rings;
	double initial_and_final_radius			   = 1.1 * libphysica::natural_units::rSun;
	unsigned int minimum_number_of_scatterings = 1;
	unsigned int maximum_number_of_scatterings = 300;
	unsigned long int maximum_free_time_steps  = 1e8;

	// Results
	unsigned long int number_of_trajectories;
	unsigned long int number_of_free_particles;
	unsigned long int number_of_reflected_particles;
	unsigned long int number_of_captured_particles;
	double average_number_of_scatterings;
	double computing_time;

	std::vector<unsigned long int> number_of_data_points;

  public:
	std::vector<std::vector<libphysica::DataPoint>> data;

	Simulation_Data(unsigned int sample_size, double u_min = 0.0, unsigned int iso_rings = 1);

	void Configure(double initial_radius, unsigned int min_scattering, unsigned int max_scattering, unsigned long int max_free_steps = 1e8);

	void Generate_Data(obscura::DM_Particle& DM, Solar_Model& solar_model, int mpi_rank = 0, unsigned int fixed_seed = 0);
	void Generate_Data_RMA(obscura::DM_Particle& DM, Solar_Model& solar_model, unsigned int fixed_seed = 0);
	void Generate_Data_Ring(obscura::DM_Particle& DM, Solar_Model& solar_model, unsigned int fixed_seed = 0);

	double Free_Ratio();
	double Capture_Ratio();
	double Reflection_Ratio();

	double Lowest_Speed(unsigned int iso_ring = 0) const;
	double Highest_Speed(unsigned int iso_ring = 0) const;

	void Print_Summary(unsigned int mpi_rank = 0);
};

#endif