#include "Reflection_Spectrum.hpp"

#include "gtest/gtest.h"
#include <mpi.h>

// Headers from libphysica
#include "Natural_Units.hpp"

// Headers from obscura
#include "DM_Distribution.hpp"
#include "DM_Particle_Standard.hpp"

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

TEST(TestReflectionSpectrum, TestSpectrum)
{
	// ARRANGE
	Solar_Model solar_model;
	double mDM = 0.1;
	obscura::DM_Particle_SI DM(mDM);
	DM.Set_Sigma_Proton(1e-1 * pb);

	solar_model.Interpolate_Total_DM_Scattering_Rate(DM, 100, 50);
	obscura::Standard_Halo_Model SHM;
	Simulation_Data data_set(10, 0, 1);
	int fixed_seed = 13;
	data_set.Generate_Data(DM, solar_model, fixed_seed);
	double v = 1e-3;
	// ACT & ASSERT
	Reflection_Spectrum spectrum(data_set, solar_model, SHM, mDM);
	EXPECT_DOUBLE_EQ(spectrum.Differential_Spectrum(v), 4.0 * M_PI * AU * AU * spectrum.Differential_DM_Flux(v));

	EXPECT_GT(spectrum.PDF_Speed(v), 0.0);

	std::function<double(double)> pdf = [&spectrum](double v) {
		return spectrum.PDF_Speed(v);
	};
	double norm = libphysica::Integrate(pdf, spectrum.Minimum_DM_Speed(), spectrum.Maximum_DM_Speed(), 1e-3);
	EXPECT_NEAR(norm, 1.0, 1e-3);

	double flux_1 = spectrum.Differential_DM_Flux(v);
	spectrum.Set_Distance(2.0 * AU);
	double flux_2 = spectrum.Differential_DM_Flux(v);
	EXPECT_DOUBLE_EQ(flux_1, 4.0 * flux_2);
}

TEST(TestReflectionSpectrum, TestDMEnteringRate)
{
	// ARRANGE
	Solar_Model solar_model;
	obscura::Standard_Halo_Model SHM;
	double mDM_1 = 1.0;
	double mDM_2 = 2.0;
	// ACT
	double rate_1 = DM_Entering_Rate(solar_model, SHM, mDM_1);
	double rate_2 = DM_Entering_Rate(solar_model, SHM, mDM_2);
	// ASSERT
	ASSERT_DOUBLE_EQ(rate_1, 2.0 * rate_2);
	ASSERT_NEAR(rate_1, 1.06689e+30 / sec, 1.0e27 / sec);
}