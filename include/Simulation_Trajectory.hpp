#ifndef __Simulation_Trajectory_hpp_
#define __Simulation_Trajectory_hpp_

#include <random>

// Headers from libphysica
#include "Natural_Units.hpp"

// Headers from obscura
#include "DM_Particle.hpp"

#include "Simulation_Utilities.hpp"
#include "Solar_Model.hpp"

// 1. Trajectory class
struct Trajectory
{
	Event current_event;
	unsigned int number_of_scatterings;
	bool continue_simulation;

	double radius_max					   = 1.5 * libphysica::natural_units::rSun;
	unsigned int number_of_scatterings_max = 300;
	unsigned long int time_steps_max	   = 1e8;

	// Constructor
	Trajectory();
	explicit Trajectory(const Event& initial_event);

	// Trajectory simulation
	void Propagate_Freely(obscura::DM_Particle& DM, Solar_Model& model, std::mt19937& PRNG);
	void Scatter(obscura::DM_Particle& DM, Solar_Model& model, std::mt19937& PRNG);
};

// 2. Numerical solution of the equations of motion (EoM)
class EoM_Solver
{
  private:
	double time, angular_momentum;
	double radius, phi, v_radial;

	libphysica::Vector axis_x, axis_y, axis_z;

	std::vector<double> tolerance = {1.0 * libphysica::natural_units::km, 1e-3 * libphysica::natural_units::km / libphysica::natural_units::sec, 1e-6};
	double dtMin				  = 1e-6 * libphysica::natural_units::sec;
	double dtMax				  = 100 * libphysica::natural_units::sec;

	double dr_dt(double v);
	double dv_dt(double r, double mass);
	double dphi_dt(double r);

  public:
	double dt;

	//Constructors
	EoM_Solver();
	explicit EoM_Solver(const Event& event);

	void Runge_Kutta_45_Step(double mass);

	double Current_Radius();
	double Current_Speed();

	bool Gravitationally_Bound(Solar_Model& model);

	Event Event_In_3D();
};

// 3. Simulation of a DM particle's orbit through the Sun
extern Trajectory Simulate_Trajectory(const Event& initial_condition, obscura::DM_Particle& DM, Solar_Model& model, std::mt19937& PRNG);

#endif