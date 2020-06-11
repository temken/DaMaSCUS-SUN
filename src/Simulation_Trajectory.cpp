#include "Simulation_Trajectory.hpp"

#include <algorithm>
#include <cmath>

// Headers from libphysica
#include "Statistics.hpp"

// Headers from obscura
#include "Astronomy.hpp"

using namespace libphysica::natural_units;

// 1. Result of one trajectory
Trajectory_Result::Trajectory_Result(const Event& event_ini, const Event& event_final, unsigned long int nScat)
: initial_event(event_ini), final_event(event_final), number_of_scatterings(nScat)
{
}

bool Trajectory_Result::Particle_Reflected() const
{
	double r	= final_event.Radius();
	double vesc = sqrt(2 * G_Newton * mSun / r);
	return r > rSun && final_event.Speed() > vesc && number_of_scatterings > 0;
}

bool Trajectory_Result::Particle_Free() const
{
	return number_of_scatterings == 0;
}

bool Trajectory_Result::Particle_Captured() const
{
	double r	= final_event.Radius();
	double vesc = sqrt(2 * G_Newton * mSun / r);
	return r < rSun || final_event.Speed() < vesc;
}

void Trajectory_Result::Print_Summary(Solar_Model& solar_model, unsigned int MPI_rank)
{
	if(MPI_rank == 0)
	{
		std::cout << SEPARATOR
				  << "Trajectory result summary" << std::endl
				  << std::endl
				  << "Number of scatterings:\t" << number_of_scatterings << std::endl
				  << "Simulation time [days]:\t" << libphysica::Round(In_Units(final_event.time, day)) << std::endl
				  << "Final radius [rSun]:\t" << libphysica::Round(In_Units(final_event.Radius(), rSun)) << std::endl
				  << "Final speed [km/sec]:\t" << libphysica::Round(In_Units(final_event.Speed(), km / sec)) << std::endl
				  << "Free particle:\t\t[" << (Particle_Free() ? "x" : " ") << "]" << std::endl
				  << "Captured:\t\t[" << (Particle_Captured() ? "x" : " ") << "]" << std::endl
				  << "Reflection:\t\t[" << (Particle_Reflected() ? "x" : " ") << "]";

		if(Particle_Reflected())
		{
			double u_i = sqrt(initial_event.Asymptotic_Speed_Sqr(solar_model));
			double u_f = sqrt(final_event.Asymptotic_Speed_Sqr(solar_model));
			std::cout << "\t(ratio u_f/u_i = " << libphysica::Round(u_f / u_i) << ")" << std::endl;
		}
		else
			std::cout << std::endl;
		std::cout << SEPARATOR << std::endl;
	}
}

// 2. Simulator
Trajectory_Simulator::Trajectory_Simulator(const Solar_Model& model)
: solar_model(model)
{
	// Pseudo-random number generator
	std::random_device rd;
	PRNG.seed(rd());
}
bool Trajectory_Simulator::Propagate_Freely(Event& current_event, obscura::DM_Particle& DM, std::ofstream& f)
{
	// 1. Define a equation-of-motion-solver in the orbital plane
	Free_Particle_Propagator particle_propagator(current_event);

	// 2. Simulate a free orbit
	double minus_log_xi			 = -log(libphysica::Sample_Uniform(PRNG));
	bool success				 = false;
	unsigned long int time_steps = 0;
	while(time_steps < maximum_time_steps && !success)
	{
		time_steps++;
		double r_before = particle_propagator.Current_Radius();
		particle_propagator.Runge_Kutta_45_Step(solar_model.Mass(r_before));
		double r_after = particle_propagator.Current_Radius();
		double v_after = particle_propagator.Current_Speed();

		if(save_trajectories && time_steps % 40 == 0)
		{
			Event event = particle_propagator.Event_In_3D();
			f << In_Units(event.time, sec) << "\t"
			  << In_Units(event.position[0], km) << "\t" << In_Units(event.position[1], km) << "\t" << In_Units(event.position[2], km) << "\t"
			  << In_Units(event.velocity[0], km / sec) << "\t" << In_Units(event.velocity[1], km / sec) << "\t" << In_Units(event.velocity[2], km / sec) << std::endl;
		}

		// Check for scatterings
		bool scattering = false;
		if(r_after < rSun)
		{
			minus_log_xi -= particle_propagator.time_step * solar_model.Total_DM_Scattering_Rate(DM, r_after, v_after);
			if(minus_log_xi < 0.0)
				scattering = true;
		}

		//Check for reflection
		bool reflection = false;
		if(r_before < maximum_distance && r_after > maximum_distance && v_after > solar_model.Local_Escape_Speed(r_after))
			reflection = true;

		if(reflection || scattering)
			success = true;
	}
	current_event = particle_propagator.Event_In_3D();
	return success;
}

int Trajectory_Simulator::Sample_Target(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(r > rSun)
	{
		std::cerr << "Error in Trajectory_Simulator::Sample_Target(): r > rSun." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		double xi		  = libphysica::Sample_Uniform(PRNG);
		double sum		  = 0.0;
		double total_rate = solar_model.Total_DM_Scattering_Rate(DM, r, DM_speed);
		//Electron
		double rate_electron = solar_model.DM_Scattering_Rate_Electron(DM, r, DM_speed);
		sum += rate_electron / total_rate;
		if(sum > xi)
			return -1;
		//Nuclei
		for(unsigned int i = 0; i < solar_model.target_isotopes.size(); i++)
		{
			double rate_nucleus = solar_model.DM_Scattering_Rate_Nucleus(DM, r, DM_speed, i);
			sum += rate_nucleus / total_rate;
			if(sum > xi)
				return i;
		}
		std::cerr << "Error in Trajectory_Simulator::Sample_Target(): No target could be sampled." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

libphysica::Vector Trajectory_Simulator::Sample_Target_Velocity(double r, double mass)
{
	if(r > rSun)
	{
		std::cerr << "Error in Trajectory_Simulator::Sample_Target_Velocity(): r > rSun." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		// 1. Sample direction
		double phi			 = libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
		double costheta		 = libphysica::Sample_Uniform(PRNG, -1.0, 1.0);
		libphysica::Vector n = libphysica::Spherical_Coordinates(1.0, acos(costheta), phi);
		// 2. Sample speed
		double temperature				  = solar_model.Temperature(r);
		double a						  = sqrt(temperature / mass);
		std::function<double(double)> cdf = [a](double v) {
			return libphysica::CDF_Maxwell_Boltzmann(v, a);
		};
		double v = libphysica::Inverse_Transform_Sampling(cdf, 0.0, 1.0, PRNG);
		return v * n;
	}
}

double Trajectory_Simulator::Sample_Scattering_Angle_Nucleus(obscura::DM_Particle& DM, Solar_Isotope& isotope)
{
	return acos(libphysica::Sample_Uniform(PRNG, -1.0, 1.0));
}

double Trajectory_Simulator::Sample_Scattering_Angle_Electron(obscura::DM_Particle& DM)
{
	return acos(libphysica::Sample_Uniform(PRNG, -1.0, 1.0));
}

libphysica::Vector Trajectory_Simulator::New_DM_Velocity(double scattering_angle, double DM_mass, double target_mass, libphysica::Vector& vel_DM, libphysica::Vector& vel_target)
{
	//Construction of n, the unit vector pointing into the direction of vfinal.
	libphysica::Vector ev = vel_DM.Normalized();
	double cosphi		  = libphysica::Sample_Uniform(PRNG, -1.0, 1.0);
	double sinphi		  = sqrt(1.0 - cosphi * cosphi);
	double cosalpha		  = cos(scattering_angle);
	double sinalpha		  = sqrt(1.0 - cosalpha * cosalpha);
	double aux			  = sqrt(1.0 - pow(ev[2], 2.0));
	libphysica::Vector n({cosalpha * ev[0] + (sinalpha * (-ev[0] * ev[2] * cosphi + ev[1] * sinphi)) / aux,
						  cosalpha * ev[1] + (sinalpha * (-ev[1] * ev[2] * cosphi - ev[0] * sinphi)) / aux,
						  cosalpha * ev[2] + aux * cosphi * sinalpha});
	double relative_speed = (vel_target - vel_DM).Norm();

	return target_mass * relative_speed / (target_mass + DM_mass) * n + (DM_mass * vel_DM + target_mass * vel_target) / (target_mass + DM_mass);
}

void Trajectory_Simulator::Scatter(Event& current_event, obscura::DM_Particle& DM)
{
	double r = current_event.Radius();
	double v = current_event.Speed();
	// 1. Find target properties.
	int target_index = Sample_Target(DM, r, v);

	double target_mass;
	if(target_index == -1)
		target_mass = mElectron;
	else
		target_mass = solar_model.target_isotopes[target_index].mass;

	libphysica::Vector vel_target = Sample_Target_Velocity(r, target_mass);

	// 2. Sample the scattering angle
	double scattering_angle;
	if(target_index == -1)
		scattering_angle = Sample_Scattering_Angle_Electron(DM);
	else
		scattering_angle = Sample_Scattering_Angle_Nucleus(DM, solar_model.target_isotopes[target_index]);

	// 3. Construct the final DM velocity
	current_event.velocity = New_DM_Velocity(scattering_angle, DM.mass, target_mass, current_event.velocity, vel_target);
}

void Trajectory_Simulator::Toggle_Trajectory_Saving(unsigned int max_trajectories)
{
	saved_trajectories	   = 0;
	saved_trajectories_max = max_trajectories;
	save_trajectories	   = !save_trajectories;
}

void Trajectory_Simulator::Fix_PRNG_Seed(int fixed_seed)
{
	PRNG.seed(fixed_seed);
}

Trajectory_Result Trajectory_Simulator::Simulate(const Event& initial_condition, obscura::DM_Particle& DM)
{
	std::ofstream f;
	if(save_trajectories && saved_trajectories < saved_trajectories_max)
	{
		saved_trajectories++;
		std::string path = TOP_LEVEL_DIR "results/";
		f.open(path + "trajectory_" + std::to_string(saved_trajectories) + ".txt");
	}
	Event current_event						= initial_condition;
	long unsigned int number_of_scatterings = 0;
	while(Propagate_Freely(current_event, DM, f) && number_of_scatterings < maximum_scatterings)
	{
		if(current_event.Radius() < rSun)
		{
			Scatter(current_event, DM);
			number_of_scatterings++;
		}
		else
			break;
	}
	if(save_trajectories)
		f.close();
	return Trajectory_Result(initial_condition, current_event, number_of_scatterings);
}

// 3. Equation of motion solution with Runge-Kutta-Fehlberg
Free_Particle_Propagator::Free_Particle_Propagator(const Event& event)
{
	// 1. Coordinate system
	axis_x = event.position.Normalized();
	axis_z = event.position.Cross(event.velocity).Normalized();
	axis_y = axis_z.Cross(axis_x);

	// 2. Coordinates
	time			 = event.time;
	radius			 = event.Radius();
	phi				 = 0.0;
	v_radial		 = (radius == 0) ? event.Speed() : event.position.Dot(event.velocity) / radius;
	angular_momentum = (event.position.Cross(event.velocity)).Dot(axis_z);
}

double Free_Particle_Propagator::dr_dt(double v)
{
	return v;
}

double Free_Particle_Propagator::dv_dt(double r, double mass)
{
	return pow(angular_momentum, 2) / pow(r, 3) - G_Newton * mass / pow(r, 2);
}

double Free_Particle_Propagator::dphi_dt(double r)
{
	return angular_momentum / pow(r, 2);
}

void Free_Particle_Propagator::Runge_Kutta_45_Step(double mass)
{
	// RK coefficients:
	double k_r[6];
	double k_v[6];
	double k_p[6];

	k_r[0] = time_step * dr_dt(v_radial);
	k_v[0] = time_step * dv_dt(radius, mass);
	k_p[0] = time_step * dphi_dt(radius);

	k_r[1] = time_step * dr_dt(v_radial + k_v[0] / 4.0);
	k_v[1] = time_step * dv_dt(radius + k_r[0] / 4.0, mass);
	// k_p[1]=	dt*dphi_dt(radius+k_r[0]/4.0,J);

	k_r[2] = time_step * dr_dt(v_radial + 3.0 / 32.0 * k_v[0] + 9.0 / 32.0 * k_v[1]);
	k_v[2] = time_step * dv_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1], mass);
	k_p[2] = time_step * dphi_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1]);

	k_r[3] = time_step * dr_dt(v_radial + 1932.0 / 2197.0 * k_v[0] - 7200.0 / 2197.0 * k_v[1] + 7296.0 / 2197.0 * k_v[2]);
	k_v[3] = time_step * dv_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2], mass);
	k_p[3] = time_step * dphi_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2]);

	k_r[4] = time_step * dr_dt(v_radial + 439.0 / 216.0 * k_v[0] - 8.0 * k_v[1] + 3680.0 / 513.0 * k_v[2] - 845.0 / 4104.0 * k_v[3]);
	k_v[4] = time_step * dv_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3], mass);
	k_p[4] = time_step * dphi_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3]);

	k_r[5] = time_step * dr_dt(v_radial - 8.0 / 27.0 * k_v[0] + 2.0 * k_v[1] - 3544.0 / 2565.0 * k_v[2] + 1859.0 / 4104.0 * k_v[3] - 11.0 / 40.0 * k_v[4]);
	k_v[5] = time_step * dv_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4], mass);
	k_p[5] = time_step * dphi_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4]);

	// New values with Runge Kutta 4 and Runge Kutta 5
	double radius_4	  = radius + 25.0 / 216.0 * k_r[0] + 1408.0 / 2565.0 * k_r[2] + 2197.0 / 4101.0 * k_r[3] - 1.0 / 5.0 * k_r[4];
	double v_radial_4 = v_radial + 25.0 / 216.0 * k_v[0] + 1408.0 / 2565.0 * k_v[2] + 2197.0 / 4101.0 * k_v[3] - 1.0 / 5.0 * k_v[4];
	double phi_4	  = phi + 25.0 / 216.0 * k_p[0] + 1408.0 / 2565.0 * k_p[2] + 2197.0 / 4101.0 * k_p[3] - 1.0 / 5.0 * k_p[4];
	double radius_5	  = radius + 16.0 / 135.0 * k_r[0] + 6656.0 / 12825.0 * k_r[2] + 28561.0 / 56430.0 * k_r[3] - 9.0 / 50.0 * k_r[4] + 2.0 / 55.0 * k_r[5];
	double v_radial_5 = v_radial + 16.0 / 135.0 * k_v[0] + 6656.0 / 12825.0 * k_v[2] + 28561.0 / 56430.0 * k_v[3] - 9.0 / 50.0 * k_v[4] + 2.0 / 55.0 * k_v[5];
	double phi_5	  = phi + 16.0 / 135.0 * k_p[0] + 6656.0 / 12825.0 * k_p[2] + 28561.0 / 56430.0 * k_p[3] - 9.0 / 50.0 * k_p[4] + 2.0 / 55.0 * k_p[5];

	// Error and adapting the time step
	std::vector<double> errors = {fabs(radius_5 - radius_4), fabs(v_radial_5 - v_radial_4), fabs(phi_5 - phi_4)};
	std::vector<double> deltas;
	for(int i = 0; i < 3; i++)
		deltas.push_back(0.84 * pow(error_tolerances[i] / errors[i], 1.0 / 4.0));
	double delta		 = *std::min_element(std::begin(deltas), std::end(deltas));
	double time_step_new = delta * time_step;
	time_step_new		 = std::max(time_step_min, std::min(time_step_new, time_step_max));

	// Check if errors fall below the tolerance
	if(errors[0] < error_tolerances[0] && errors[1] < error_tolerances[1] && errors[2] < error_tolerances[2])
	{
		time	  = time + time_step;
		radius	  = radius_4;
		v_radial  = v_radial_4;
		phi		  = phi_4;
		time_step = time_step_new;
	}
	else
	{
		time_step = time_step_new;
		Runge_Kutta_45_Step(mass);
	}
}

double Free_Particle_Propagator::Current_Time()
{
	return time;
}

double Free_Particle_Propagator::Current_Radius()
{
	return radius;
}

double Free_Particle_Propagator::Current_Speed()
{
	if(radius == 0 || angular_momentum == 0)
		return v_radial;
	else
		return sqrt(v_radial * v_radial + angular_momentum * angular_momentum / radius / radius);
}

Event Free_Particle_Propagator::Event_In_3D()
{
	double v_phi			= angular_momentum / pow(radius, 2);
	libphysica::Vector xNew = radius * (cos(phi) * axis_x + sin(phi) * axis_y);
	libphysica::Vector vNew = (v_radial * cos(phi) - v_phi * radius * sin(phi)) * axis_x + (v_radial * sin(phi) + radius * v_phi * cos(phi)) * axis_y;

	return Event(time, xNew, vNew);
}