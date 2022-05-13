#include "gtest/gtest.h"

#include "Medium_Effects.hpp"

using namespace DaMaSCUS_SUN;
using namespace std::complex_literals;

TEST(TestMediumEffects, TestPlasmaDispersionFunction)
{
	// ARRANGE
	double x					= 0.0;
	std::complex<double> result = 0.0 + 1.0i * std::sqrt(M_PI);
	// ACT & ASSERT
	std::complex<double> z = Plasma_Dispersion_Function(x);
	EXPECT_NEAR(z.real(), result.real(), 1e-6);
	EXPECT_NEAR(z.imag(), result.imag(), 1e-6);
}