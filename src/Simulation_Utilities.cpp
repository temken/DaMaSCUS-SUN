#include "Simulation_Utilities.hpp"

#include <cmath>

#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/DM_Halo_Models.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

// 1. Event Class
Event::Event()
: time(0.0), position({0, 0, 0}), velocity({0, 0, 0})
{
}

Event::Event(double t, const libphysica::Vector& pos, const libphysica::Vector& vel)
: time(t), position(pos), velocity(vel)
{
}

double Event::Radius() const
{
	return position.Norm();
}

double Event::Speed() const
{
	return velocity.Norm();
}

double Event::Angular_Momentum() const
{
	return position.Cross(velocity).Norm();
}

double Event::Asymptotic_Speed_Sqr(Solar_Model& solar_model) const
{
	double r	 = Radius();
	double v	 = Speed();
	double vesc	 = solar_model.Local_Escape_Speed(r);
	double u_sqr = v * v - vesc * vesc;
	return u_sqr;
}

double Event::Isoreflection_Angle(const libphysica::Vector& vel_sun) const
{
	return acos(position.Normalized().Dot(vel_sun.Normalized()));
}

int Event::Isoreflection_Ring(const libphysica::Vector& vel_sun, unsigned int number_of_rings) const
{
	double theta					= Isoreflection_Angle(vel_sun);
	std::vector<double> ring_angles = Isoreflection_Ring_Angles(number_of_rings);
	for(unsigned int ring = 0; ring < ring_angles.size(); ring++)
		if(theta <= ring_angles[ring])
			return ring;
	std::cerr << "Error in Event::Isoreflection_Ring(): Angle = " << theta << " out of bound." << std::endl;
	std::exit(EXIT_FAILURE);
}

Event Event::In_Units(double unit_distance, double unit_time) const
{
	return Event(libphysica::natural_units::In_Units(time, unit_time), libphysica::natural_units::In_Units(position, unit_distance), libphysica::natural_units::In_Units(velocity, unit_distance / unit_time));
}

// Overload <<
std::ostream& operator<<(std::ostream& output, const Event& event)
{
	return output << "{"
				  << event.time
				  << ","
				  << event.position
				  << ","
				  << event.velocity
				  << "}";
}

// 2. Generator of initial conditions
double PDF_Initial_Speed(double v, obscura::DM_Distribution& halo_model, Solar_Model& solar_model)
{
	double v_esc			 = solar_model.Local_Escape_Speed(rSun);
	double v_average		 = halo_model.Average_Speed();
	double v_inverse_average = halo_model.Eta_Function(0.0);
	return halo_model.PDF_Speed(v) * (v + v_esc * v_esc / v) / (v_average + v_esc * v_esc * v_inverse_average);
}

// Conditional pdf for cos_theta given a speed value v
double PDF_Cos_Theta(double cos_theta, double v, obscura::DM_Distribution& halo_model)
{
	double normalization = 2.0 * M_PI * v * v / halo_model.PDF_Speed(v);

	// Construct velocity vector with angle theta to the Sun's velocity
	libphysica::Vector vel_sun = dynamic_cast<obscura::Standard_Halo_Model*>(&halo_model)->Get_Observer_Velocity();
	libphysica::Vector vel	   = libphysica::Spherical_Coordinates(v, acos(cos_theta), 0.0, vel_sun);

	return normalization * halo_model.PDF_Velocity(vel);
}

Event Initial_Conditions(obscura::DM_Distribution& halo_model, Solar_Model& solar_model, std::mt19937& PRNG)
{
	// 1. Initial velocity
	// 1.1. Sample initial speed u asymptotically far from the Sun.
	std::function<double(double)> pdf_v = [&halo_model, &solar_model](double v) {
		return PDF_Initial_Speed(v, halo_model, solar_model);
	};
	double u = libphysica::Rejection_Sampling(pdf_v, halo_model.Minimum_DM_Speed(), halo_model.Maximum_DM_Speed(), 1000.0, PRNG);

	// 1.2. Sample cos(theta) where theta is the angle between v and v_sun.
	std::function<double(double)> pdf_cos_theta = [u, &halo_model](double cos_theta) {
		return PDF_Cos_Theta(cos_theta, u, halo_model);
	};
	double y_max	 = PDF_Cos_Theta(-1.0, u, halo_model);
	double cos_theta = libphysica::Rejection_Sampling(pdf_cos_theta, -1.0, 1.0, y_max, PRNG);

	// 1.3. Construct velocity vector
	libphysica::Vector vel_sun			= dynamic_cast<obscura::Standard_Halo_Model*>(&halo_model)->Get_Observer_Velocity();
	double phi							= libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	libphysica::Vector initial_velocity = libphysica::Spherical_Coordinates(u, acos(cos_theta), phi, vel_sun);

	// 1.4. Blue-shift the speed
	double asymptotic_distance = 1000.0 * AU;
	double vesc_asymptotic	   = solar_model.Local_Escape_Speed(asymptotic_distance);
	double v				   = sqrt(u * u + vesc_asymptotic * vesc_asymptotic);
	initial_velocity		   = v * initial_velocity.Normalized();

	// 2. Initial position
	// 2.1 Find the maximum impact parameter such that the particle still hits the Sun.
	double v_esc				= solar_model.Local_Escape_Speed(rSun);
	double impact_parameter_max = sqrt(u * u + v_esc * v_esc) / v * rSun;
	libphysica::Vector e_z		= (-1.0) * initial_velocity.Normalized();
	libphysica::Vector e_x({0, e_z[2], -e_z[1]});
	e_x.Normalize();
	libphysica::Vector e_y = e_z.Cross(e_x);

	// 2.2 Find a random point in the plane.
	double phi_disk						= libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	double xi							= libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
	double impact_parameter				= sqrt(xi) * impact_parameter_max;
	libphysica::Vector initial_position = asymptotic_distance * e_z + impact_parameter * (cos(phi_disk) * e_x + sin(phi_disk) * e_y);

	return Event(0.0, initial_position, initial_velocity);
}

// 3. Analytically propagate a particle at event on a hyperbolic Kepler orbit to a radius R (without passing the periapsis)
void Hyperbolic_Kepler_Shift(Event& event, double R_final)
{
	// 1. Initial event
	double R_initial		= event.Radius();
	double v_initial		= event.Speed();
	double angular_momentum = event.Angular_Momentum();

	if(R_final < rSun || R_initial < rSun)
	{
		std::cerr << "Error in Hyperbolic_Kepler_Shift(): Orbits inside the Sun cannot be described analytically." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// 2. Asymptotic speed
	double vEsc = sqrt(2 * G_Newton * mSun / R_initial);
	double u2	= v_initial * v_initial - vEsc * vEsc;

	// 3. Kepler orbit parameter
	double semi_major_axis	= G_Newton * mSun / u2;
	double semilatus_rectum = angular_momentum * angular_momentum / G_Newton / mSun;
	double eccentricity		= sqrt(1.0 + semilatus_rectum / semi_major_axis);
	// double perihelion		= semi_major_axis * ( eccentricity - 1.0);

	// 4. Initial and final orbital angle
	double theta_initial = libphysica::Sign(R_final - R_initial) * acos(1.0 / eccentricity * (semilatus_rectum / R_initial - 1.0));
	double theta_final	 = libphysica::Sign(R_final - R_initial) * acos(1.0 / eccentricity * (semilatus_rectum / R_final - 1.0));

	// 5. Axis vectors of orbital coordinate system
	libphysica::Vector axis_z = event.position.Cross(event.velocity).Normalized();
	libphysica::Vector axis_x = cos(theta_initial) * event.position.Normalized() + sin(theta_initial) * event.position.Normalized().Cross(axis_z);
	libphysica::Vector axis_y = axis_z.Cross(axis_x);

	// 6. Final time, position, and velocity
	// 6.1 Time
	// double F1 = acosh((eccentricity + cos(theta_initial)) / (1.0 + eccentricity * cos(theta_initial)));
	// double M1 = eccentricity * sinh(F1) - F1;
	// double t1 = sqrt(pow(+semi_major_axis, 3) / G_Newton / mSun) * M1;
	// double F2 = acosh((eccentricity + cos(theta_final)) / (1.0 + eccentricity * cos(theta_final)));
	// double M2 = eccentricity * sinh(F2) - F2;
	// double t2 = sqrt(pow(+semi_major_axis, 3) / G_Newton / mSun) * M2;
	// event.time += libphysica::Sign(R_final - R_initial) * (t2 - t1);
	// 6.2 Position and Velocity
	event.position = R_final * cos(theta_final) * axis_x + R_final * sin(theta_final) * axis_y;
	event.velocity = sqrt(G_Newton * mSun / semilatus_rectum) * (eccentricity * sin(theta_final) * event.position.Normalized() + (1.0 + eccentricity * cos(theta_final)) * axis_z.Cross(event.position.Normalized()));
}

// 4. Equiareal isodetection rings
std::vector<double> Isoreflection_Ring_Angles(unsigned int number_of_rings)
{
	std::vector<double> thetas;
	double theta = 0;
	for(unsigned int i = 0; i < number_of_rings - 1; i++)
	{
		theta = acos(cos(theta) - 2.0 / number_of_rings);
		thetas.push_back(theta);
	}
	thetas.push_back(180. * deg);
	return thetas;
}

}	// namespace DaMaSCUS_SUN