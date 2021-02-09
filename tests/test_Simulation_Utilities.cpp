#include "gtest/gtest.h"

// Headers from obscura
#include "Astronomy.hpp"

#include "Simulation_Trajectory.hpp"
#include "Simulation_Utilities.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

// 1. Event class
TEST(TestSimulationUtilities, TestEventConstructor)
{
	// ARRANGE
	Event default_event;
	double t = 7 * sec;
	libphysica::Vector x({1 * meter, 2 * meter, 3 * meter});
	libphysica::Vector v({4 * meter / sec, 5 * meter / sec, 6 * meter / sec});
	Event event(t, x, v);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(default_event.time, 0.0);
	for(int i = 0; i < 3; i++)
		EXPECT_DOUBLE_EQ(default_event.position[i], libphysica::Vector({0, 0, 0})[i]);
	for(int i = 0; i < 3; i++)
		EXPECT_DOUBLE_EQ(default_event.velocity[i], libphysica::Vector({0, 0, 0})[i]);
	EXPECT_DOUBLE_EQ(In_Units(event.time, sec), 7.0);
	for(int i = 0; i < 3; i++)
		EXPECT_DOUBLE_EQ(event.position[i], x[i]);
	for(int i = 0; i < 3; i++)
		EXPECT_DOUBLE_EQ(event.velocity[i], v[i]);
}

TEST(TestSimulationUtilities, TestRadius)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector x({1, 2, 3});
	libphysica::Vector v({4, 5, 6});
	Event event(t, x, v);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(event.Radius(), sqrt(14));
}

TEST(TestSimulationUtilities, TestSpeed)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector x({1, 2, 3});
	libphysica::Vector v({4, 5, 6});
	Event event(t, x, v);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(event.Speed(), sqrt(77));
}

TEST(TestSimulationUtilities, TestAngularMomentum)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector x({1, 2, 3});
	libphysica::Vector v({4, 5, 6});
	Event event(t, x, v);
	double J = x.Cross(v).Norm();
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(event.Angular_Momentum(), J);
}

TEST(TestSimulationUtilities, TestAsymptoticSpeedSqr)
{
	// ARRANGE
	Solar_Model SSM;
	double r	= 5.0 * rSun;
	double u	= 100.0 * km / sec;
	double vesc = SSM.Local_Escape_Speed(r);
	double v	= sqrt(u * u + vesc * vesc);

	libphysica::Vector x({0, r, 0});
	libphysica::Vector vel({v, 0, 0});
	Event event(0, x, vel);

	double tol = 1.0e-10;
	// ACT & ASSERT
	ASSERT_NEAR(event.Asymptotic_Speed_Sqr(SSM), u * u, tol);
}

TEST(TestSimulationUtilities, TestIsoreflectionAngle)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector x({1, 2, 3});
	libphysica::Vector v({4, 5, 6});
	Event event(t, x, v);

	libphysica::Vector vSun = obscura::Sun_Velocity();
	double phi				= acos(x.Normalized() * vSun.Normalized());
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(event.Isoreflection_Angle(vSun), phi);
}

TEST(TestSimulationUtilities, TestIsoreflectionRing)
{
	// ARRANGE
	libphysica::Vector vSun = obscura::Sun_Velocity();
	unsigned int rings		= 3;
	double t				= 0;
	libphysica::Vector x({1, 2, 3});
	libphysica::Vector v({4, 5, 6});
	Event event(t, x, v);

	// ACT & ASSERT
	event.position = vSun;
	EXPECT_EQ(event.Isoreflection_Ring(vSun, rings), 0);
	event.position = -1.0 * vSun;
	EXPECT_EQ(event.Isoreflection_Ring(vSun, rings), 2);
	event.position = libphysica::Vector({vSun[1], -vSun[0], 0});
	EXPECT_EQ(event.Isoreflection_Ring(vSun, rings), 1);
}

TEST(TestSimulationUtilities, TestInUnits)
{
	// ARRANGE
	double t = 7 * sec;
	libphysica::Vector x({1 * km, 2 * km, 3 * km});
	libphysica::Vector v({4 * km / sec, 5 * km / sec, 6 * km / sec});
	Event event(t, x, v);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(event.In_Units(km, minute).time, 7.0 / 60);
	EXPECT_DOUBLE_EQ(event.In_Units(rSun, sec).position[0], km / rSun);
	EXPECT_DOUBLE_EQ(event.In_Units(km, sec).position[1], 2);
	EXPECT_DOUBLE_EQ(event.In_Units(km, sec).position[2], 3);
	EXPECT_DOUBLE_EQ(event.In_Units(km, sec).velocity[0], 4);
	EXPECT_DOUBLE_EQ(event.In_Units(km, sec).velocity[1], 5);
	EXPECT_DOUBLE_EQ(event.In_Units(km, sec).velocity[2], 6);
}

// 2. Generator of initial conditions

TEST(TestSimulationUtilities, TestInitialConditions)
{
	// ARRANGE
	std::random_device rd;
	std::mt19937 PRNG;
	PRNG.seed(3);

	Solar_Model SSM;

	obscura::Standard_Halo_Model SHM;
	SHM.Set_Observer_Velocity(libphysica::Vector({0, 0, 0}));

	double R_distance	= 500 * AU;
	unsigned int trials = 100;
	// ACT & ASSERT
	for(unsigned int i = 0; i < trials; i++)
	{
		Event IC	= Initial_Conditions(SHM, SSM, PRNG);
		double r	= IC.Radius();
		double v	= IC.Speed();
		double vesc = SSM.Local_Escape_Speed(r);
		double Jmax = rSun * sqrt(v * v + SSM.Local_Escape_Speed(rSun) * SSM.Local_Escape_Speed(rSun));
		ASSERT_GE(IC.Radius(), R_distance);
		ASSERT_GE(IC.Speed(), vesc);
		ASSERT_LE(IC.Angular_Momentum(), Jmax);
		ASSERT_GE(IC.Angular_Momentum(), 0.0);
	}
}

// 3. Analytically propagate a particle at event on a hyperbolic Kepler orbit to a radius R (without passing the periapsis)
TEST(TestSimulationUtilities, TestHyperbolicKeplerShift)
{
	// ARRANGE
	// std::random_device rd;
	// std::mt19937 PRNG(rd());
	int fixed_seed = 123;
	std::mt19937 PRNG(fixed_seed);

	Solar_Model SSM;

	obscura::Standard_Halo_Model SHM;
	SHM.Set_Observer_Velocity(libphysica::Vector({0, 0, 0}));

	int trials = 100;
	for(int k = 0; k < trials; k++)
	{
		// ARRANGE
		Event IC = Initial_Conditions(SHM, SSM, PRNG);
		Hyperbolic_Kepler_Shift(IC, 5.0 * rSun);

		Free_Particle_Propagator eom(IC);
		while(eom.Current_Radius() > rSun)
			eom.Runge_Kutta_45_Step(mSun);
		Event x_ref = eom.Event_In_3D();
		// ACT
		Hyperbolic_Kepler_Shift(IC, rSun);
		//ASSERT
		for(int i = 0; i < 3; i++)
			ASSERT_NEAR(IC.position[i], x_ref.position[i], 0.001 * rSun);
		for(int i = 0; i < 3; i++)
			ASSERT_NEAR(IC.velocity[i], x_ref.velocity[i], km / sec);
		ASSERT_LE(IC.Speed(), x_ref.Speed());
	}
}

// 4. Equiareal isodetection rings
TEST(TestSimulationUtilities, TestIsoreflectionRingAngles)
{
	// ARRANGE
	std::vector<double> two_rings	  = {90 * deg, 180 * deg};
	std::vector<double> hundred_rings = {36.8699 * deg, 53.1301 * deg, 66.4218 * deg, 78.463 * deg, 90. * deg, 101.537 * deg, 113.578 * deg, 126.87 * deg, 143.13 * deg, 180. * deg};

	double tol = 1.0e-3 * deg;
	// ACT & ASSERT
	for(unsigned int i = 0; i < two_rings.size(); i++)
		EXPECT_NEAR(Isoreflection_Ring_Angles(2)[i], two_rings[i], tol);
	for(unsigned int i = 0; i < hundred_rings.size(); i++)
		EXPECT_NEAR(Isoreflection_Ring_Angles(10)[i], hundred_rings[i], tol);
}