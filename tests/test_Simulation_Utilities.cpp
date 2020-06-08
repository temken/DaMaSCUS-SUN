#include "gtest/gtest.h"

// Headers from obscura
#include "Astronomy.hpp"

#include "Simulation_Trajectory.hpp"
#include "Simulation_Utilities.hpp"

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

TEST(TestSimulationUtilities, TestInKmSec)
{
	// ARRANGE
	double t = 7 * sec;
	libphysica::Vector x({1 * km, 2 * km, 3 * km});
	libphysica::Vector v({4 * km / sec, 5 * km / sec, 6 * km / sec});
	Event event(t, x, v);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(event.In_km_sec().time, 7);
	EXPECT_DOUBLE_EQ(event.In_km_sec().position[0], 1);
	EXPECT_DOUBLE_EQ(event.In_km_sec().position[1], 2);
	EXPECT_DOUBLE_EQ(event.In_km_sec().position[2], 3);
	EXPECT_DOUBLE_EQ(event.In_km_sec().velocity[0], 4);
	EXPECT_DOUBLE_EQ(event.In_km_sec().velocity[1], 5);
	EXPECT_DOUBLE_EQ(event.In_km_sec().velocity[2], 6);
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
	std::random_device rd;
	std::mt19937 PRNG(rd());
	// PRNG.seed(3);
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;
	SHM.Set_Observer_Velocity(libphysica::Vector({0, 0, 0}));

	int trials = 500;
	for(int k = 0; k < trials; k++)
	{
		// ARRANGE
		Event IC = Initial_Conditions(SHM, SSM, PRNG);
		Hyperbolic_Kepler_Shift(IC, 5.0 * rSun);

		EoM_Solver eom(IC);
		while(eom.Current_Radius() > rSun)
			eom.Runge_Kutta_45_Step(mSun);
		Event x_ref = eom.Event_In_3D();
		// ACT
		Hyperbolic_Kepler_Shift(IC, rSun);
		//ASSERT
		for(int i = 0; i < 3; i++)
			ASSERT_NEAR(IC.position[i], x_ref.position[i], 0.005 * rSun);
		for(int i = 0; i < 3; i++)
			ASSERT_NEAR(IC.velocity[i], x_ref.velocity[i], km / sec);
		ASSERT_LE(IC.Speed(), x_ref.Speed());
	}
}

// 4. Equiareal isodetection rings
