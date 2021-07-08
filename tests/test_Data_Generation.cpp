#include "Data_Generation.hpp"

#include "gtest/gtest.h"
#include <mpi.h>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

int main(int argc, char* argv[])
{
	int result = 0;

	::testing::InitGoogleTest(&argc, argv);
	MPI_Init(&argc, &argv);
	result = RUN_ALL_TESTS();
	MPI_Finalize();
	return result;
}

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
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	unsigned int sample_size = 2;

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 100, 50);

	// ACT
	Simulation_Data data_set(sample_size);
	data_set.Generate_Data(DM, SSM, SHM);

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
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 100, 50);

	unsigned int sample_size = 10;
	// ACT
	Simulation_Data data_set(sample_size);
	data_set.Configure(1.1 * rSun, 0, 500);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT
	ASSERT_DOUBLE_EQ(data_set.Free_Ratio(), 1.0);
}

TEST(TestDataGeneration, TestDataSetCaptureRatio)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(1.0 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 100, 50);

	unsigned int sample_size = 50;

	// ACT
	Simulation_Data data_set(sample_size);
	data_set.Configure(1.1 * rSun, 0, 500);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT
	ASSERT_GE(data_set.Capture_Ratio(), 0.0);
}

TEST(TestDataGeneration, TestDataSetReflectionRatio)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(1.0 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 100, 50);

	unsigned int sample_size = 10;

	// ACT
	Simulation_Data data_set(sample_size);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT
	ASSERT_GT(data_set.Reflection_Ratio(), 0.0);
}

TEST(TestDataGeneration, TestSpeedFunctions)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(1.0 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 100, 50);

	unsigned int sample_size = 10;
	double u_min			 = 0.0001;
	// ACT
	Simulation_Data data_set(sample_size, u_min);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT
	EXPECT_DOUBLE_EQ(data_set.Minimum_Speed(), 0.75 * u_min);
	EXPECT_GT(data_set.Lowest_Speed(), data_set.Minimum_Speed());
	EXPECT_GT(data_set.Highest_Speed(), data_set.Lowest_Speed());
}

// TEST(TestDataGeneration, TestDataSetPrintSummary)
// {
// 	// ARRANGE
// 	// void Print_Summary(unsigned int mpi_rank = 0);
// 	// ACT & ASSERT
// }
