#include "gtest/gtest.h"

#include "Medium_Effects.hpp"

using namespace DaMaSCUS_SUN;
using namespace std::complex_literals;

TEST(TestMediumEffects, TestDawsonIntegral)
{
	// ARRANGE
	std::vector<double> xs		= {0.1, 0.5, 2.3, 10.7, 22.0};
	std::vector<double> results = {0.099336, 0.424436, 0.249053, 0.0469358, 0.0227508};
	// ACT & ASSERT
	for(unsigned int i = 0; i < xs.size(); i++)
		EXPECT_NEAR(Dawson_Integral(xs[i]), results[i], 1e-6);
}

TEST(TestMediumEffects, TestErfi)
{
	// ARRANGE
	std::vector<double> xs		= {0.15, 0.45, 2.35, 10.24, 22.9};
	std::vector<double> results = {0.170535, 0.544232, 68.3326, 1.91557e44, 1.38157e226};
	// ACT & ASSERT
	for(unsigned int i = 0; i < xs.size(); i++)
		EXPECT_NEAR(Erfi(xs[i]), results[i], 1e-5 * results[i]);
}

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