#ifndef __Solar_Model_hpp_
#define __Solar_Model_hpp_

#include <random>

// Headers from libphysica
#include "Linear_Algebra.hpp"
#include "Numerics.hpp"

// Headers from obscura
#include "Target_Nucleus.hpp"

//1. Nuclear targets in the Sun
	class Solar_Nucleus : public obscura::Element
	{
	private:
		libphysica::Interpolation number_density;
	public:
		Solar_Nucleus(obscura::Isotope isotope, std::vector<std::vector<double>>& density_table);
		Solar_Nucleus(std::vector<obscura::Isotope> isotopes, std::vector<std::vector<double>>& density_table);

		double Number_Density(double r);

		libphysica::Vector Sample_Velocity(double temperature, std::mt19937& PRNG);
	};

//2. Solar model
	class Solar_Model
	{
	private:
		libphysica::Interpolation mass, temperature, escape_velocity_squared;

	public:
		Solar_Model();

		double Temperature(double r);
		double Escape_Velocity(double r);
		double Mass(double r);

		void Print_Summary(int MPI_rank = 0);
	};

#endif
