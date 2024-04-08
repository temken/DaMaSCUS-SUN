#ifndef __Simulation_Trajectory_hpp_
#define __Simulation_Trajectory_hpp_

#include <fstream>
#include <random>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Particle.hpp"

#include "Simulation_Utilities.hpp"
#include "Solar_Model.hpp"

namespace DaMaSCUS_SUN
{

// 1. Result of one trajectory
struct Trajectory_Result
{
	Event initial_event, final_event;
	unsigned long int number_of_scatterings;
	double radius_last_scattering;
	double radius_deepest_scattering;

	Trajectory_Result(const Event& event_ini, const Event& event_final, unsigned long int nScat, double r_last = -1.0, double r_deepest = -1.0);

	bool Particle_Reflected() const;
	bool Particle_Free() const;
	bool Particle_Captured(Solar_Model& solar_model) const;

	void Print_Summary(Solar_Model& solar_model, unsigned int mpi_rank = 0);
};

// 2. Simulator
class Trajectory_Simulator
{
  private:
	Solar_Model solar_model;

	unsigned int saved_trajectories, saved_trajectories_max;
	bool save_trajectories = false;
	double v_max		   = 0.75;

	bool Propagate_Freely(Event& current_event, obscura::DM_Particle& DM, std::ofstream& f);

	int Sample_Target(obscura::DM_Particle& DM, double r, double vDM);
	libphysica::Vector Sample_Momentum_Transfer(int target_index, obscura::DM_Particle& DM, const libphysica::Vector& DM_velocity, double r);

  public:
	std::mt19937 PRNG;
	unsigned long int maximum_time_steps;
	unsigned int maximum_scatterings;
	double maximum_distance;

	Trajectory_Simulator(const Solar_Model& model, unsigned long int max_time_steps = 1e8, unsigned int max_scatterings = 500, double max_distance = 1.1 * libphysica::natural_units::rSun);

	void Toggle_Trajectory_Saving(unsigned int max_trajectories = 50);
	void Fix_PRNG_Seed(int fixed_seed);

	void Scatter(Event& current_event, obscura::DM_Particle& DM);
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

  public:
	double time_step = 0.1 * libphysica::natural_units::sec;

	explicit Free_Particle_Propagator(const Event& event);

	void Runge_Kutta_45_Step(double mass);

	double Current_Time();
	double Current_Radius();
	double Current_Speed();

	Event Event_In_3D();
};

}	// namespace DaMaSCUS_SUN

#endif