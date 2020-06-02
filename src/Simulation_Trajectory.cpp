#include "Simulation_Trajectory.hpp"

#include <algorithm>
#include <cmath>

// Headers from obscura
#include "Astronomy.hpp"

using namespace libphysica::natural_units;

// 1. Trajectory class
Trajectory::Trajectory()
: current_event(Event()), number_of_scatterings(0), continue_simulation(true)
{
}

Trajectory::Trajectory(const Event& initial_event)
: current_event(initial_event), number_of_scatterings(0), continue_simulation(true)
{
}

void Trajectory::Propagate_Freely(obscura::DM_Particle& DM, Solar_Model& model, std::mt19937& PRNG)
{
	// 1. Define a equation-of-motion-solver in the orbital plane
	EoM_Solver eom(current_event);

	// 2. Simulate a free orbit
	unsigned long int time_steps = 0;
	while(continue_simulation)
	{
		time_steps++;
		double r_before = eom.Current_Radius();
		eom.Runge_Kutta_45_Step(model.Mass(r_before));
		double r_after = eom.Current_Radius();

		bool solar_reflection = r_before < radius_max && r_after > radius_max && eom.Gravitationally_Bound(model) == false;
		if(solar_reflection || time_steps > time_steps_max)
			continue_simulation = false;
	}
	current_event = eom.Event_In_3D();
}

void Trajectory::Scatter(obscura::DM_Particle& DM, Solar_Model& model, std::mt19937& PRNG)
{
	number_of_scatterings++;
	if(number_of_scatterings > number_of_scatterings_max)
		continue_simulation = false;
	else
	{
	}
}

// 2. Numerical solution of the equations of motion (EoM)
double EoM_Solver::dr_dt(double v)
{
	return v;
}

double EoM_Solver::dv_dt(double r, double mass)
{
	return pow(angular_momentum, 2) / pow(r, 3) - G_Newton * mass / pow(r, 2);
}

double EoM_Solver::dphi_dt(double r)
{
	return angular_momentum / pow(r, 2);
}

EoM_Solver::EoM_Solver()
: time(0), angular_momentum(0), radius(0), phi(0), v_radial(0), axis_x({1, 0, 0}), axis_y({0, 1, 0}), axis_z({0, 0, 1}), dt(0.1 * sec)
{
}

EoM_Solver::EoM_Solver(const Event& event)
: time(event.time), radius(event.Radius()), phi(0), dt(0.1 * sec)
{
	axis_x = event.position.Normalized();
	axis_z = event.position.Cross(event.velocity).Normalized();
	axis_y = axis_z.Cross(axis_x);

	angular_momentum = (event.position.Cross(event.velocity)).Dot(axis_z);
	v_radial		 = event.position.Dot(event.velocity) / radius;
}

void EoM_Solver::Runge_Kutta_45_Step(double mass)
{
	//RK coefficients:
	double k_r[6];
	double k_v[6];
	double k_p[6];

	k_r[0] = dt * dr_dt(v_radial);
	k_v[0] = dt * dv_dt(radius, mass);
	k_p[0] = dt * dphi_dt(radius);

	k_r[1] = dt * dr_dt(v_radial + k_v[0] / 4.0);
	k_v[1] = dt * dv_dt(radius + k_r[0] / 4.0, mass);
	// k_p[1]=	dt*dphi_dt(radius+k_r[0]/4.0,J);

	k_r[2] = dt * dr_dt(v_radial + 3.0 / 32.0 * k_v[0] + 9.0 / 32.0 * k_v[1]);
	k_v[2] = dt * dv_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1], mass);
	k_p[2] = dt * dphi_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1]);

	k_r[3] = dt * dr_dt(v_radial + 1932.0 / 2197.0 * k_v[0] - 7200.0 / 2197.0 * k_v[1] + 7296.0 / 2197.0 * k_v[2]);
	k_v[3] = dt * dv_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2], mass);
	k_p[3] = dt * dphi_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2]);

	k_r[4] = dt * dr_dt(v_radial + 439.0 / 216.0 * k_v[0] - 8.0 * k_v[1] + 3680.0 / 513.0 * k_v[2] - 845.0 / 4104.0 * k_v[3]);
	k_v[4] = dt * dv_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3], mass);
	k_p[4] = dt * dphi_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3]);

	k_r[5] = dt * dr_dt(v_radial - 8.0 / 27.0 * k_v[0] + 2.0 * k_v[1] - 3544.0 / 2565.0 * k_v[2] + 1859.0 / 4104.0 * k_v[3] - 11.0 / 40.0 * k_v[4]);
	k_v[5] = dt * dv_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4], mass);
	k_p[5] = dt * dphi_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4]);

	//Compute new values and errors
	double r4				= radius + 25.0 / 216.0 * k_r[0] + 1408.0 / 2565.0 * k_r[2] + 2197.0 / 4101.0 * k_r[3] - 1.0 / 5.0 * k_r[4];
	double v4				= v_radial + 25.0 / 216.0 * k_v[0] + 1408.0 / 2565.0 * k_v[2] + 2197.0 / 4101.0 * k_v[3] - 1.0 / 5.0 * k_v[4];
	double phi4				= phi + 25.0 / 216.0 * k_p[0] + 1408.0 / 2565.0 * k_p[2] + 2197.0 / 4101.0 * k_p[3] - 1.0 / 5.0 * k_p[4];
	double r5				= radius + 16.0 / 135.0 * k_r[0] + 6656.0 / 12825.0 * k_r[2] + 28561.0 / 56430.0 * k_r[3] - 9.0 / 50.0 * k_r[4] + 2.0 / 55.0 * k_r[5];
	double v5				= v_radial + 16.0 / 135.0 * k_v[0] + 6656.0 / 12825.0 * k_v[2] + 28561.0 / 56430.0 * k_v[3] - 9.0 / 50.0 * k_v[4] + 2.0 / 55.0 * k_v[5];
	double phi5				= phi + 16.0 / 135.0 * k_p[0] + 6656.0 / 12825.0 * k_p[2] + 28561.0 / 56430.0 * k_p[3] - 9.0 / 50.0 * k_p[4] + 2.0 / 55.0 * k_p[5];
	std::vector<double> err = {fabs(r5 - r4), fabs(v5 - v4), fabs(phi5 - phi4)};

	//New stepsize
	std::vector<double> deltas;
	for(int i = 0; i < 3; i++)
		deltas.push_back(0.84 * pow(tolerance[i] / err[i], 1.0 / 4.0));
	double delta = *std::min_element(std::begin(deltas), std::end(deltas));
	//Next steps
	if(err[0] < tolerance[0] && err[1] < tolerance[1] && err[2] < tolerance[2])
	{
		time += dt;
		radius	 = r4;
		v_radial = v4;
		phi		 = phi4;
		dt		 = std::max(dtMin, std::min(delta * dt, dtMax));
	}
	else
	{
		dt = std::max(dtMin, std::min(delta * dt, dtMax));
		Runge_Kutta_45_Step(mass);
	}
}

double EoM_Solver::Current_Radius()
{
	return radius;
}

double EoM_Solver::Current_Speed()
{
	return sqrt(v_radial * v_radial + angular_momentum * angular_momentum / radius / radius);
}

bool EoM_Solver::Gravitationally_Bound(Solar_Model& model)
{
	return Current_Speed() < model.Local_Escape_Speed(radius);
}

Event EoM_Solver::Event_In_3D()
{
	double v_phi			= angular_momentum / pow(radius, 2);
	libphysica::Vector xNew = radius * (cos(phi) * axis_x + sin(phi) * axis_y);
	libphysica::Vector vNew = (v_radial * cos(phi) - v_phi * radius * sin(phi)) * axis_x + (v_radial * sin(phi) + radius * v_phi * cos(phi)) * axis_y;

	return Event(time, xNew, vNew);
}

// 3. Simulation of a DM particle's orbit through the Sun
Trajectory Simulate_Trajectory(const Event& initial_condition, obscura::DM_Particle& DM, Solar_Model& model, std::mt19937& PRNG)
{
	Trajectory traj(initial_condition);
	while(traj.continue_simulation)
	{
		traj.Propagate_Freely(DM, model, PRNG);
		if(traj.continue_simulation)
			traj.Scatter(DM, model, PRNG);
	}
	return traj;
}
