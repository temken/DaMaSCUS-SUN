#ifndef __Simulation_Trajectory_hpp_
#define __Simulation_Trajectory_hpp_

#include <fstream>
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
	Event initial_event, final_event;
	unsigned long int number_of_scatterings;

	Trajectory_Result(const Event& event_ini, const Event& event_final, unsigned long int nScat);

	bool Particle_Reflected() const;
	bool Particle_Free() const;
	bool Particle_Captured() const;

	void Print_Summary(Solar_Model& solar_model, unsigned int MPI_rank = 0);
};

// 2. Simulator
class Trajectory_Simulator
{
  private:
	Solar_Model solar_model;

	unsigned int saved_trajectories, saved_trajectories_max;
	bool save_trajectories = false;

	bool Propagate_Freely(Event& current_event, obscura::DM_Particle& DM, std::ofstream& f);

	int Sample_Target(obscura::DM_Particle& DM, double r, double DM_speed);
	libphysica::Vector Sample_Target_Velocity(double r, double mass);
	double Sample_Scattering_Angle_Nucleus(obscura::DM_Particle& DM, Solar_Isotope& isotope);
	double Sample_Scattering_Angle_Electron(obscura::DM_Particle& DM);
	libphysica::Vector New_DM_Velocity(double scattering_angle, double DM_mass, double target_mass, libphysica::Vector& vel_DM, libphysica::Vector& vel_target);

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
	std::vector<double> error_tolerances = {1.0 * libphysica::natural_units::km, 1.0e-3 * libphysica::natural_units::km / libphysica::natural_units::sec, 1.0e-7};
	double time_step_min				 = 1.0e-6 * libphysica::natural_units::sec;
	double time_step_max				 = 1.0e2 * libphysica::natural_units::sec;

  public:
	double time_step = 0.1 * libphysica::natural_units::sec;

	explicit Free_Particle_Propagator(const Event& event);

	void Runge_Kutta_45_Step(double mass);

	double Current_Time();
	double Current_Radius();
	double Current_Speed();

	Event Event_In_3D();
};

#endif