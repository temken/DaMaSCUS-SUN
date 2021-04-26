#include "gtest/gtest.h"

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Particle_Standard.hpp"

#include "Simulation_Trajectory.hpp"
#include "Simulation_Utilities.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

// 1. Result of one trajectory
TEST(TestSimulationTrajectory, TestTrajectoryResultConstructor)
{
	// ARRANGE
	double t_i = 0;
	libphysica::Vector r_i({1, 2, 3});
	libphysica::Vector v_i({4, 5, 6});
	double t_f = 7;
	libphysica::Vector r_f({8, 9, 10});
	libphysica::Vector v_f({11, 12, 13});
	Event initial_event(t_i, r_i, v_i);
	Event final_event(t_f, r_f, v_f);
	unsigned int nscat = 14;
	// ACT
	Trajectory_Result result(initial_event, final_event, nscat);
	// ASSERT
	EXPECT_EQ(result.initial_event.time, t_i);
	EXPECT_EQ(result.initial_event.position, r_i);
	EXPECT_EQ(result.initial_event.velocity, v_i);
	EXPECT_EQ(result.final_event.time, t_f);
	EXPECT_EQ(result.final_event.position, r_f);
	EXPECT_EQ(result.final_event.velocity, v_f);
	EXPECT_EQ(result.number_of_scatterings, nscat);
}

TEST(TestSimulationTrajectory, TestParticleReflected)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector r_in({0, 0.5 * rSun, 0});
	libphysica::Vector r_out({0, 1.5 * rSun, 0});
	libphysica::Vector v_1({0, 0, 0});
	libphysica::Vector v_2({0, 1000 * km / sec, 0});
	Event initial_event;
	Event not_reflected_1(t, r_in, v_1);
	Event not_reflected_2(t, r_in, v_2);
	Event not_reflected_3(t, r_out, v_1);
	Event reflected_1(t, r_out, v_2);

	// ACT & ASSERT
	EXPECT_FALSE(Trajectory_Result(initial_event, not_reflected_1, 1).Particle_Reflected());
	EXPECT_FALSE(Trajectory_Result(initial_event, not_reflected_2, 1).Particle_Reflected());
	EXPECT_FALSE(Trajectory_Result(initial_event, not_reflected_3, 1).Particle_Reflected());
	EXPECT_TRUE(Trajectory_Result(initial_event, reflected_1, 1).Particle_Reflected());
	EXPECT_FALSE(Trajectory_Result(initial_event, reflected_1, 0).Particle_Reflected());
}

TEST(TestSimulationTrajectory, TestParticleFree)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector r({0, 0.5 * rSun, 0});
	libphysica::Vector v({0, 1000 * km / sec, 0});
	Event initial_event;
	Event final_event(t, r, v);

	// ACT & ASSERT
	EXPECT_FALSE(Trajectory_Result(initial_event, final_event, 1).Particle_Free());
	EXPECT_TRUE(Trajectory_Result(initial_event, final_event, 0).Particle_Free());
}

TEST(TestSimulationTrajectory, TestParticleCaptured)
{
	// ARRANGE
	double t = 0;
	libphysica::Vector r_in({0, 0.5 * rSun, 0});
	libphysica::Vector r_out({0, 1.5 * rSun, 0});
	libphysica::Vector v_1({0, 0, 0});
	libphysica::Vector v_2({0, 1000 * km / sec, 0});
	Event initial_event;
	Event captured_1(t, r_in, v_1);
	Event captured_2(t, r_in, v_2);
	Event captured_3(t, r_out, v_1);
	Event not_captured(t, r_out, v_2);

	// ACT & ASSERT
	EXPECT_TRUE(Trajectory_Result(initial_event, captured_1, 0).Particle_Captured());
	EXPECT_TRUE(Trajectory_Result(initial_event, captured_2, 0).Particle_Captured());
	EXPECT_TRUE(Trajectory_Result(initial_event, captured_3, 0).Particle_Captured());
	EXPECT_FALSE(Trajectory_Result(initial_event, not_captured, 0).Particle_Captured());
}

// TEST(TestSimulationTrajectory, TestResultPrintSummary)
// {
// 	//ARRANGE
// 	//ACT
// 	//ASSERT
// }

// 2. Simulator

TEST(TestSimulationTrajectory, TestScatter)
{
	// ARRANGE
	obscura::DM_Particle_SI DM(0.5 * GeV);
	DM.Set_Sigma_Proton(pb);
	double time = 0.0;
	libphysica::Vector r({0.25 * rSun, 0.25 * rSun, 0.25 * rSun});
	libphysica::Vector vel({100 * km / sec, 100 * km / sec, 100 * km / sec});
	Event event(time, r, vel);
	Solar_Model SSM;
	Trajectory_Simulator simulator(SSM);
	// ACT & ASSERT
	int trials = 10;
	for(int i = 0; i < trials; i++)
	{
		libphysica::Vector v_ini = event.velocity;
		simulator.Scatter(event, DM);
		for(int j = 0; j < 3; j++)
			EXPECT_NE(event.velocity[j], v_ini[j]);
	}
}

TEST(TestSimulationTrajectory, TestSimulate)
{
	// ARRANGE
	obscura::DM_Particle_SI DM(0.5 * GeV);
	DM.Set_Sigma_Proton(0.1 * pb);
	Solar_Model SSM;
	Trajectory_Simulator simulator(SSM);
	obscura::Standard_Halo_Model SHM;
	// ACT & ASSERT
	int trials = 5;
	for(int i = 0; i < trials; i++)
	{
		Event IC = Initial_Conditions(SHM, SSM, simulator.PRNG);
		Hyperbolic_Kepler_Shift(IC, 1.5 * rSun);
		Trajectory_Result result = simulator.Simulate(IC, DM);
		if(result.Particle_Reflected() || result.Particle_Free())
			ASSERT_NEAR(result.final_event.Radius(), simulator.maximum_distance, 0.001 * rSun);
		else
			ASSERT_GT(result.number_of_scatterings, 0);
	}
}

// TEST(TestSimulationTrajectory, TestSimulatorPrintSummary)
// {
// 	//ARRANGE
// 	//ACT
// 	//ASSERT
// }

// 3. Equation of motion solution with Runge-Kutta-Fehlberg

TEST(TestSimulationTrajectory, TestRungeKuttaStep)
{
	// ASSERT
	double t = 0;
	libphysica::Vector r({0.25 * rSun, 0.5 * rSun, 0.25 * rSun});
	libphysica::Vector v({km / sec, 1000 * km / sec, km / sec});
	Event event(t, r, v);
	Free_Particle_Propagator propagator(event);
	double initial_timestep = propagator.time_step;
	// ACT
	propagator.Runge_Kutta_45_Step(mSun);
	// // ASSERT
	EXPECT_NE(propagator.time_step, initial_timestep);
	EXPECT_GT(propagator.Current_Time(), event.time);
	EXPECT_GT(propagator.Current_Radius(), event.Radius());
	EXPECT_LT(propagator.Current_Speed(), event.Speed());
}

TEST(TestSimulationTrajectory, TestPropagatorCurrentRadius)
{
	// ASSERT
	double t	  = 13 * sec;
	double radius = 17 * km;
	libphysica::Vector r({radius, 0, 0});
	libphysica::Vector vel;
	Event event(t, r, vel);
	Free_Particle_Propagator propagator(event);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(propagator.Current_Radius(), radius);
}

TEST(TestSimulationTrajectory, TestPropagatorCurrentSpeed)
{
	// ASSERT
	double t = 13 * sec;
	double v = 123.0 * km / sec;
	libphysica::Vector r;
	libphysica::Vector vel({0, v, 0});
	Event event(t, r, vel);
	Free_Particle_Propagator propagator(event);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(propagator.Current_Speed(), v);
}

TEST(TestSimulationTrajectory, TestPropagatorEventIn3D)
{
	// ASSERT
	double t = 13 * sec;
	libphysica::Vector r({km, 0.5 * rSun, km});
	libphysica::Vector v({km / sec, 1000 * km / sec, km / sec});
	Event event(t, r, v);
	Free_Particle_Propagator propagator(event);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(propagator.Event_In_3D().time, t);
	for(int i = 0; i < 3; i++)
		EXPECT_DOUBLE_EQ(propagator.Event_In_3D().position[i], r[i]);
	for(int i = 0; i < 3; i++)
		EXPECT_DOUBLE_EQ(propagator.Event_In_3D().velocity[i], v[i]);
}
