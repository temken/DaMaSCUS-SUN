#ifndef __Data_Generation_hpp_
#define __Data_Generation_hpp_

#include <vector>

// Headers from libphysica
#include "Statistics.hpp"

// Headers from obscura
#include "DM_Particle.hpp"

#include "Simulation_Trajectory.hpp"

struct Simulation_Data_Set
{
	unsigned long int number_of_trajectories;
	unsigned long int number_of_free_particles;
	unsigned long int number_of_reflected_particles;
	unsigned long int number_of_captured_particles;
	double simulation_time;

	std::vector<std::vector<libphysica::DataPoint>> data;

	Simulation_Data_Set();

	unsigned int Sampe_Size(unsigned int isoreflection_ring = 0);

	double Free_Ratio();
	double Capture_Ratio();
	double Reflection_Ratio();

	void Print_Summary(unsigned int MPI_rank = 0);
};

extern Simulation_Data_Set Generate_Data(unsigned int sample_size, obscura::DM_Particle& DM, Solar_Model& solar_model, double u_min = 0.0, unsigned int isoreflection_rings = 1);

#endif