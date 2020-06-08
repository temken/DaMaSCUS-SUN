#ifndef __Simulation_Utilities_hpp_
#define __Simulation_Utilities_hpp_

#include <random>

// Headers from libphysica
#include "Linear_Algebra.hpp"
#include "Natural_Units.hpp"

// Headers from obscura
#include "DM_Distribution.hpp"

#include "Solar_Model.hpp"

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

	double Isoreflection_Angle(const libphysica::Vector& vel_sun) const;
	int Isoreflection_Ring(const libphysica::Vector& vel_sun, unsigned int number_of_rings) const;

	Event In_km_sec() const;

	//Overloading the output operator <<
	friend std::ostream& operator<<(std::ostream& output, const Event& event);
};

// 2. Generator of initial conditions
extern Event Initial_Conditions(obscura::DM_Distribution& halo_model, Solar_Model& model, std::mt19937& PRNG);

// 3. Analytically propagate a particle at event on a hyperbolic Kepler orbit to a radius R (without passing the periapsis)
extern void Hyperbolic_Kepler_Shift(Event& event, double R_final);

// 4. Equiareal isodetection rings
extern std::vector<double> Isoreflection_Ring_Angles(unsigned int number_of_rings);

#endif