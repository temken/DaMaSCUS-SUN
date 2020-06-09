#ifndef __Simulation_Trajectory_hpp_
#define __Simulation_Trajectory_hpp_

#include <functional>
#include <random>

// Headers from libphysica
#include "Natural_Units.hpp"

// Headers from obscura
#include "DM_Particle.hpp"

#include "Simulation_Utilities.hpp"
#include "Solar_Model.hpp"

// 1. Result of one trajectory
struct Trajectory_Result
{
	Event final_event;
	unsigned long int number_of_scatterings;

	Trajectory_Result(const Event& event, unsigned long int nScat);

	bool Particle_Reflected() const;
	bool Particle_Free() const;
	bool Particle_Captured() const;

	void Print_Summary(unsigned int MPI_rank = 0);
};

// 2. Simulator
class Trajectory_Simulator
{
  private:
	std::mt19937 PRNG;
	Solar_Model solar_model;

	bool Propagate_Freely(Event& current_event, obscura::DM_Particle& DM);

	void Scatter(Event& current_event, obscura::DM_Particle& DM);

  public:
	unsigned long int maximum_time_steps  = 1e8;
	unsigned long int maximum_scatterings = 300;
	double maximum_distance				  = 1.5 * libphysica::natural_units::rSun;
	bool save_trajectory_to_file		  = false;

	Trajectory_Simulator(const Solar_Model& model);

	Trajectory_Result Simulate(const Event& initial_condition, obscura::DM_Particle& DM);
};

// 3. Equation of motion solution with Runge-Kutta-Fehlberg
class Free_Particle_Propagator
{
  private:
	double time, radius, phi, v_radial;
	double angular_momentum;
	libphysica::Vector axis_x, axis_y, axis_z;

	double dr_dt(double v);
	double dv_dt(double r, double mass);
	double dphi_dt(double r);
	std::vector<double> error_tolerances;
	double time_step_min, time_step_max;
	double time_step;

  public:
	explicit Free_Particle_Propagator(const Event& event);

	void Runge_Kutta_45_Step(double mass);

	double Current_Radius();
	double Current_Speed();

	Event Event_In_3D();
};

#endif