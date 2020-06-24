#include "Data_Generation.hpp"

#include "gtest/gtest.h"

// Headers from libphysica
#include "Natural_Units.hpp"

// Headers from obscura
#include "DM_Particle_Standard.hpp"

using namespace libphysica::natural_units;

TEST(TestDataGeneration, TestDataSetConstructor)
{
	// ARRANGE
	unsigned int sample_size = 100;
	double u_min			 = 0.0;
	unsigned int iso_rings	 = 10;
	// ACT
	Simulation_Data data_set(sample_size, u_min, iso_rings);
	// ASSERT
	ASSERT_EQ(data_set.data.size(), iso_rings);
	for(auto& set : data_set.data)
		ASSERT_EQ(set.size(), 0);
}

TEST(TestDataGeneration, TestGenerateData)
{
	// ARRANGE
	Solar_Model SSM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	unsigned int sample_size = 2;

	// ACT
	Simulation_Data data_set(sample_size);
	data_set.Generate_Data(DM, SSM);

	// ASSERT
	ASSERT_EQ(data_set.data[0].size(), sample_size);
}

TEST(TestDataGeneration, TestConfigure)
{
	// ARRANGE
	unsigned int sample_size		 = 2;
	double r						 = 2.0 * rSun;
	unsigned int min_scattering		 = 1;
	unsigned int max_scattering		 = 1;
	unsigned long int max_time_steps = 1e4;
	// ACT
	Simulation_Data data_set(sample_size);
	data_set.Configure(r, min_scattering, max_scattering, max_time_steps);

	// ASSERT
	// ASSERT_EQ(data_set.data[0].size(), sample_size);
}

TEST(TestDataGeneration, TestDataFreeRatio)
{
	// ARRANGE
	Solar_Model SSM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	unsigned int sample_size = 10;

	// ACT
	Simulation_Data data_set(sample_size);
	data_set.Configure(1.1 * rSun, 0, 500);
	data_set.Generate_Data(DM, SSM);

	// ASSERT
	ASSERT_DOUBLE_EQ(data_set.Free_Ratio(), 1.0);
}

TEST(TestDataGeneration, TestDataSetCaptureRatio)
{
	// ARRANGE
	Solar_Model SSM;

	obscura::DM_Particle_SI DM(1.0 * GeV);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	unsigned int sample_size = 20;

	// ACT
	Simulation_Data data_set(sample_size);
	data_set.Configure(1.1 * rSun, 0, 500);
	data_set.Generate_Data(DM, SSM);

	// ASSERT
	ASSERT_GT(data_set.Capture_Ratio(), 0.0);
}

TEST(TestDataGeneration, TestDataSetReflectionRatio)
{
	// ARRANGE
	Solar_Model SSM;

	obscura::DM_Particle_SI DM(1.0 * GeV);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	unsigned int sample_size = 20;

	// ACT
	Simulation_Data data_set(sample_size);
	data_set.Generate_Data(DM, SSM);

	// ASSERT
	ASSERT_GT(data_set.Reflection_Ratio(), 0.0);
}

// TEST(TestDataGeneration, TestDataSetPrintSummary)
// {
// 	// ARRANGE
// 	// void Print_Summary(unsigned int MPI_rank = 0);
// 	// ACT & ASSERT
// }
