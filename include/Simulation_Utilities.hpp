#ifndef __Simulation_Utilities_hpp_
#define __Simulation_Utilities_hpp_

#include <random>

#include "libphysica/Linear_Algebra.hpp"
#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Distribution.hpp"

#include "Solar_Model.hpp"

namespace DaMaSCUS_SUN
{

// 1. Event class
struct Event
{
	double time;
	libphysica::Vector position, velocity;

	// Constructors
	Event();
	Event(double t, const libphysica::Vector& pos, const libphysica::Vector& vel);

	double Radius() const;
	double Speed() const;
	double Angular_Momentum() const;
	double Asymptotic_Speed_Sqr(Solar_Model& solar_model) const;

	double Isoreflection_Angle(const libphysica::Vector& vel_sun) const;
	int Isoreflection_Ring(const libphysica::Vector& vel_sun, unsigned int number_of_rings) const;

	Event In_Units(double unit_distance, double unit_time) const;

	//Overloading the output operator <<
	friend std::ostream& operator<<(std::ostream& output, const Event& event);
};

// 2. Generator of initial conditions
extern Event Initial_Conditions(obscura::DM_Distribution& halo_model, Solar_Model& model, std::mt19937& PRNG);

// 3. Analytically propagate a particle at event on a hyperbolic Kepler orbit to a radius R (without passing the periapsis)
extern void Hyperbolic_Kepler_Shift(Event& event, double R_final);

// 4. Equiareal isoreflection rings
extern std::vector<double> Isoreflection_Ring_Angles(unsigned int number_of_rings);

}	// namespace DaMaSCUS_SUN

#endif