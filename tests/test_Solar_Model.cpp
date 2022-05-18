#include "gtest/gtest.h"

#include <cmath>
#include <mpi.h>
#include <random>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/DM_Particle_Standard.hpp"

#include "Solar_Model.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;
using namespace std::complex_literals;

int main(int argc, char* argv[])
{
	int result = 0;

	::testing::InitGoogleTest(&argc, argv);
	MPI_Init(&argc, &argv);
	result = RUN_ALL_TESTS();
	MPI_Finalize();
	return result;
}

TEST(TestSolarModel, TestPlasmaDispersionFunction)
{
	// ARRANGE
	double x					= 0.0;
	std::complex<double> result = 0.0 + 1.0i * std::sqrt(M_PI);
	// ACT & ASSERT
	std::complex<double> z = Plasma_Dispersion_Function(x);
	EXPECT_NEAR(z.real(), result.real(), 1e-6);
	EXPECT_NEAR(z.imag(), result.imag(), 1e-6);
}

TEST(TestSolarModel, TestSolarNucleus)
{
	// ARRANGE
	double rho	= 1.505e+02 * gram / cm / cm / cm;
	double f_H1 = 0.36209;
	double n_H1 = f_H1 * rho / mProton;
	Solar_Model SSM;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(SSM.target_isotopes[0].Number_Density(0.0), n_H1);
	ASSERT_DOUBLE_EQ(SSM.target_isotopes[0].Number_Density(0.001 * rSun), n_H1);
	ASSERT_DOUBLE_EQ(SSM.target_isotopes[0].Number_Density(0.00150 * rSun), n_H1);
}

TEST(TestSolarModel, TestConstructor)
{
	// ARRANGE
	// ACT
	Solar_Model SSM;
	// ASSERT
	ASSERT_EQ(SSM.name, "Standard Solar Model AGSS09");
	ASSERT_NE(SSM.target_isotopes.size(), 0);
}

TEST(TestSolarModel, TestMass)
{
	// ARRANGE
	Solar_Model SSM;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(SSM.Mass(0), 0.0);
	ASSERT_DOUBLE_EQ(SSM.Mass(rSun), mSun);
	ASSERT_DOUBLE_EQ(SSM.Mass(0.00150 * rSun), 0.0000004 * mSun);
}

TEST(TestSolarModel, TestTemperature)
{
	// ARRANGE
	Solar_Model SSM;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(SSM.Temperature(0), 1.549e+07 * Kelvin);
	EXPECT_DOUBLE_EQ(SSM.Temperature(rSun), 5800.0 * Kelvin);
	EXPECT_DOUBLE_EQ(SSM.Temperature(0.96950 * rSun), 1.544e+05 * Kelvin);
}

TEST(TestSolarModel, TestLocalEscapeSpeed)
{
	// ARRANGE
	Solar_Model SSM;
	// ACT & ASSERT
	ASSERT_NEAR(SSM.Local_Escape_Speed(0), 1384.13 * km / sec, 0.1 * km / sec);
	ASSERT_DOUBLE_EQ(SSM.Local_Escape_Speed(rSun), sqrt(2.0 * G_Newton * mSun / rSun));
	ASSERT_DOUBLE_EQ(SSM.Local_Escape_Speed(2.0 * rSun), sqrt(2.0 * G_Newton * mSun / 2.0 / rSun));
}

TEST(TestSolarModel, TestNumberDensityNucleus)
{
	// ARRANGE
	Solar_Model SSM;
	double r1 = 0.0 * rSun;
	double r2 = 0.5 * rSun;
	double r3 = 1.5 * rSun;
	// ACT & ASSERT
	for(unsigned int i = 0; i < SSM.target_isotopes.size(); i++)
	{
		EXPECT_DOUBLE_EQ(SSM.Number_Density_Nucleus(r1, i), SSM.target_isotopes[i].Number_Density(r1));
		EXPECT_DOUBLE_EQ(SSM.Number_Density_Nucleus(r2, i), SSM.target_isotopes[i].Number_Density(r2));
		EXPECT_DOUBLE_EQ(SSM.Number_Density_Nucleus(r3, i), 0.0);
	}
}

TEST(TestSolarModel, TestNumberDensityElectron)
{
	// ARRANGE
	Solar_Model SSM;
	double r		 = 0.5 * rSun;
	double nElectron = 0.0;
	for(auto& isotope : SSM.all_isotopes)
		nElectron += isotope.Z * isotope.Number_Density(r);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(SSM.Number_Density_Electron(r), nElectron);
	EXPECT_DOUBLE_EQ(SSM.Number_Density_Electron(1.5 * rSun), 0.0);
}

TEST(TestSolarModel, TestDMScatteringRateElectron)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::DM_Particle_SI DM;
	DM.Set_Sigma_Electron(pb);
	double v_DM = 1e-3;
	double r1	= 0.5 * rSun;
	double r2	= 1.5 * rSun;
	// ACT & ASSERT
	EXPECT_GE(SSM.DM_Scattering_Rate_Electron(DM, r1, v_DM), 0.0);
	EXPECT_DOUBLE_EQ(SSM.DM_Scattering_Rate_Electron(DM, r2, v_DM), 0.0);
	DM.Set_Sigma_Electron(0.0);
	EXPECT_DOUBLE_EQ(SSM.DM_Scattering_Rate_Electron(DM, r1, v_DM), 0.0);
}

TEST(TestSolarModel, TestDMScatteringRateNucleus)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::DM_Particle_SI DM;
	DM.Set_Sigma_Proton(pb);
	double v_DM = 1e-3;
	double r1	= 0.5 * rSun;
	double r2	= 1.5 * rSun;
	// ACT & ASSERT
	EXPECT_GE(SSM.DM_Scattering_Rate_Nucleus(DM, r1, v_DM, 0), 0.0);
	EXPECT_DOUBLE_EQ(SSM.DM_Scattering_Rate_Nucleus(DM, r2, v_DM, 0), 0.0);
	DM.Unfix_Coupling_Ratios();
	DM.Set_Sigma_Proton(0.0);
	EXPECT_DOUBLE_EQ(SSM.DM_Scattering_Rate_Nucleus(DM, r1, v_DM, 0), 0.0);
}

TEST(TestSolarModel, TestTotalDMScatteringRate)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::DM_Particle_SI DM;
	DM.Set_Sigma_Proton(pb);
	double v_DM	 = 1e-3;
	double r1	 = 0.5 * rSun;
	double r2	 = 1.5 * rSun;
	double total = SSM.DM_Scattering_Rate_Electron(DM, r1, v_DM);
	for(unsigned int i = 0; i < SSM.target_isotopes.size(); i++)
		total += SSM.DM_Scattering_Rate_Nucleus(DM, r1, v_DM, i);
	// ACT & ASSERT
	EXPECT_GE(SSM.Total_DM_Scattering_Rate(DM, r1, v_DM), 0.0);
	EXPECT_DOUBLE_EQ(SSM.Total_DM_Scattering_Rate(DM, r2, v_DM), 0.0);
	EXPECT_DOUBLE_EQ(SSM.Total_DM_Scattering_Rate(DM, r1, v_DM), total);
	DM.Unfix_Coupling_Ratios();
	DM.Set_Sigma_Proton(0.0);
	DM.Set_Sigma_Neutron(0.0);
	DM.Set_Sigma_Electron(0.0);
	EXPECT_DOUBLE_EQ(SSM.Total_DM_Scattering_Rate(DM, r1, v_DM), 0.0);
}

TEST(TestSolarModel, TestTotalDMScatteringRateInterpolation)
{
	// ARRANGE
	int fixed_seed = 998;
	std::mt19937 PRNG(fixed_seed);
	Solar_Model SSM;
	obscura::DM_Particle_SI DM(0.01);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(pb);
	int trials		 = 500;
	double tolerance = 0.1;
	// ACT
	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 1000, 1000);
	// ASSERT
	for(int i = 0; i < trials; i++)
	{
		double r			 = libphysica::Sample_Uniform(PRNG, 0, rSun);
		double w			 = libphysica::Sample_Uniform(PRNG, 0, 0.3);
		double correct_value = SSM.Total_DM_Scattering_Rate_Computed(DM, r, w);
		double max_deviation = tolerance * correct_value;
		ASSERT_NEAR(SSM.Total_DM_Scattering_Rate(DM, r, w), correct_value, max_deviation);
	}
}

TEST(TestSolarModel, TestPrintSummary)
{
	// ARRANGE
	Solar_Model SSM;

	// ACT & ASSERT
	SSM.Print_Summary();
}
