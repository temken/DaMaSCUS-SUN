#include "gtest/gtest.h"

#include <cmath>

// Headers from libphysica
#include "Natural_Units.hpp"

#include "Solar_Model.hpp"

using namespace libphysica::natural_units;

TEST(TestSolarModel, TestSolarNucleus)
{
	// ARRANGE
	double rho	= 1.505e+02 * gram / cm / cm / cm;
	double f_H1 = 0.36209;
	double n_H1 = f_H1 * rho / mProton;
	Solar_Model SSM;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(SSM.nuclear_targets[0].Number_Density(0.0), n_H1);
	ASSERT_DOUBLE_EQ(SSM.nuclear_targets[0].Number_Density(0.001 * rSun), n_H1);
	ASSERT_DOUBLE_EQ(SSM.nuclear_targets[0].Number_Density(0.00150 * rSun), n_H1);
}

TEST(TestSolarModel, TestConstructor)
{
	// ARRANGE
	// ACT
	Solar_Model SSM;
	// ASSERT
	ASSERT_EQ(SSM.name, "Standard Solar Model AGSS09");
	ASSERT_NE(SSM.nuclear_targets.size(), 0);
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
	ASSERT_DOUBLE_EQ(SSM.Temperature(0), 1.549e+07 * Kelvin);
	ASSERT_DOUBLE_EQ(SSM.Temperature(rSun), 7.063e+04 * Kelvin);
	ASSERT_DOUBLE_EQ(SSM.Temperature(0.96950 * rSun), 1.544e+05 * Kelvin);
}

TEST(TestSolarModel, TestLocalEscapeSpeed)
{
	// ARRANGE
	Solar_Model SSM;
	// ACT & ASSERT
	ASSERT_NEAR(SSM.Local_Escape_Speed(0), 1384.6 * km / sec, 0.1 * km / sec);
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
	for(unsigned int i = 0; i < SSM.nuclear_targets.size(); i++)
	{
		EXPECT_DOUBLE_EQ(SSM.Number_Density_Nucleus(r1, i), SSM.nuclear_targets[i].Number_Density(r1));
		EXPECT_DOUBLE_EQ(SSM.Number_Density_Nucleus(r2, i), SSM.nuclear_targets[i].Number_Density(r2));
		EXPECT_DOUBLE_EQ(SSM.Number_Density_Nucleus(r3, i), 0.0);
	}
}

TEST(TestSolarModel, TestNumberDensityElectron)
{
	// ARRANGE
	Solar_Model SSM;
	double r		 = 0.5 * rSun;
	double nElectron = 0.0;
	for(auto& nucleus : SSM.nuclear_targets)
	{
		nElectron += nucleus[0].Z * nucleus.Number_Density(r);
	}
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(SSM.Number_Density_Electron(r), nElectron);
	EXPECT_DOUBLE_EQ(SSM.Number_Density_Electron(1.5 * rSun), 0.0);
}

// TEST(TestSolarModel, TestPrintSummary)
// {
// 	// ARRANGE
// 	// ACT & ASSERT
// }
