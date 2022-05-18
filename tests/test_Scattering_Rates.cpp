#include "gtest/gtest.h"

#include "libphysica/Natural_Units.hpp"

#include "Scattering_Rates.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

// 3. Thermal average of relative speed between a particle of speed v_DM and a solar thermal target.
TEST(TestSolarModel, TestThermalAveragedRelativeSpeed)
{
	// ARRANGE
	double temperature = 1.0e7 * Kelvin;
	double mass_target = mProton;
	double v_DM		   = 1e-3;
	// ACT & ASSERT
	EXPECT_NEAR(Thermal_Averaged_Relative_Speed(0.01 * Kelvin, mass_target, v_DM), v_DM, 1e-10);
	EXPECT_NEAR(Thermal_Averaged_Relative_Speed(temperature, mass_target, v_DM), 0.0017928, 1.0e-7);
}
