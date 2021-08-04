#include "gtest/gtest.h"

#include "Dark_Photon.hpp"

#include <random>

#include "libphysica/Natural_Units.hpp"

#include "obscura/Target_Nucleus.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

TEST(TestDarkPhoton, TestConstructor1)
{
	// ARRANGE
	DM_Particle_Dark_Photon DM;
	//ACT & ASSERT
	EXPECT_DOUBLE_EQ(DM.mass, 10.0);
	EXPECT_DOUBLE_EQ(DM.Sigma_Proton(), 1.0e-40 * cm * cm);
}

TEST(TestDarkPhoton, TestConstructor2)
{
	// ARRANGE
	double mDM = 100 * MeV;
	DM_Particle_Dark_Photon DM(mDM);
	//ACT & ASSERT
	EXPECT_DOUBLE_EQ(DM.mass, mDM);
	EXPECT_DOUBLE_EQ(DM.Sigma_Proton(), 1e-40 * cm * cm);
}

TEST(TestDarkPhoton, TestConstructor3)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double sigma_p = pb;
	DM_Particle_Dark_Photon DM(mDM, sigma_p);
	//ACT & ASSERT
	EXPECT_DOUBLE_EQ(DM.mass, mDM);
	EXPECT_DOUBLE_EQ(DM.Sigma_Proton(), sigma_p);
}

TEST(TestDarkPhoton, TestSigmaRatio)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double sigma_p = pb;
	DM_Particle_Dark_Photon DM(mDM, sigma_p);
	double sigma_e = pow(libphysica::Reduced_Mass(mDM, mElectron) / libphysica::Reduced_Mass(mDM, mProton), 2.0) * sigma_p;
	//ACT & ASSERT
	ASSERT_DOUBLE_EQ(DM.Sigma_Electron(), sigma_e);
}

TEST(TestDarkPhoton, TestSetMass)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double sigma_p = pb;
	DM_Particle_Dark_Photon DM(mDM, sigma_p);
	double sigma_e = DM.Sigma_Electron();
	// ACT
	DM.Set_Mass(10 * MeV);
	// ASSERT
	EXPECT_DOUBLE_EQ(DM.mass, 10 * MeV);
	EXPECT_DOUBLE_EQ(DM.Sigma_Electron(), sigma_e);
}

TEST(TestDarkPhoton, TestEpsilon)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double sigma_p = pb;
	DM_Particle_Dark_Photon DM(mDM, sigma_p);
	// ACT
	double epsilon = DM.Get_Epsilon();
	DM.Set_Epsilon(2.0 * epsilon);
	// ASSERT
	ASSERT_DOUBLE_EQ(DM.Sigma_Proton(), 4.0 * sigma_p);
}

TEST(TestDarkPhoton, TestInteractionParameter)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double sigma_p = pb;
	DM_Particle_Dark_Photon DM(mDM, sigma_p);
	//ACT & ASSERT
	double sig = 10 * pb;
	DM.Set_Interaction_Parameter(sig, "Nuclei");
	EXPECT_DOUBLE_EQ(DM.Get_Interaction_Parameter("Nuclei"), sig);
	EXPECT_DOUBLE_EQ(DM.Sigma_Proton(), sig);
	DM.Set_Interaction_Parameter(sig, "Electrons");
	EXPECT_DOUBLE_EQ(DM.Get_Interaction_Parameter("Electrons"), sig);
	EXPECT_DOUBLE_EQ(DM.Sigma_Electron(), sig);
}

TEST(TestDarkPhoton, TestSetSigma)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double sigma_p = pb;
	DM_Particle_Dark_Photon DM(mDM, sigma_p);
	//ACT & ASSERT
	double sig = 10 * pb;
	DM.Set_Sigma_Proton(sig);
	EXPECT_DOUBLE_EQ(DM.Sigma_Proton(), sig);
	DM.Set_Sigma_Electron(sig);
	EXPECT_DOUBLE_EQ(DM.Sigma_Electron(), sig);
}

TEST(TestDarkPhoton, TestFormFactors)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double sigma_p = pb;
	DM_Particle_Dark_Photon DM(mDM, sigma_p);
	double qref = aEM * mElectron;
	double q	= keV;
	double vDM	= 1e-3;
	double mMed = GeV;
	//ACT & ASSERT
	double dodq2 = DM.dSigma_dq2_Electron(q, vDM);
	DM.Set_FormFactor_DM("Contact");
	EXPECT_DOUBLE_EQ(DM.dSigma_dq2_Electron(q, vDM), dodq2);
	DM.Set_FormFactor_DM("General", mMed);
	EXPECT_DOUBLE_EQ(DM.dSigma_dq2_Electron(q, vDM), pow((qref * qref + mMed * mMed) / (q * q + mMed * mMed), 2.0) * dodq2);
	DM.Set_FormFactor_DM("Long-Range");
	EXPECT_DOUBLE_EQ(DM.dSigma_dq2_Electron(q, vDM), pow((qref * qref + mMed * mMed) / qref / qref, 2.0) * pow(qref / q, 4.0) * dodq2);
}

TEST(TestDarkPhoton, TestSigmaTotalGeneralNucleus)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double vDM	   = 1e-3;
	double r	   = 0.5 * rSun;
	double sigma_p = pb;
	DM_Particle_Dark_Photon DM(mDM, sigma_p);
	auto target = obscura::Get_Isotope(2, 4);
	// ACT
	double sigma_contact = DM.Sigma_Total_Nucleus(target, vDM, r);
	DM.Set_FormFactor_DM("General", 1e6);
	DM.Set_Sigma_Proton(sigma_p);
	double sigma_general = DM.Sigma_Total_Nucleus(target, vDM, r);
	double tol			 = 1e-6 * sigma_general;
	// ASSERT
	EXPECT_NEAR(sigma_contact, sigma_general, tol);
}

TEST(TestDarkPhoton, TestSigmaTotalGeneralElectron)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double vDM	   = 1e-3;
	double r	   = 0.5 * rSun;
	double sigma_e = pb;
	DM_Particle_Dark_Photon DM(mDM, pb);
	DM.Set_Sigma_Electron(sigma_e);
	// ACT
	double sigma_contact = DM.Sigma_Total_Electron(vDM, r);
	DM.Set_FormFactor_DM("General", 1e6);
	DM.Set_Sigma_Electron(sigma_e);
	double sigma_general = DM.Sigma_Total_Electron(vDM, r);
	double tol			 = 1e-6 * sigma_general;
	// ASSERT
	EXPECT_NEAR(sigma_contact, sigma_general, tol);
}

TEST(TestDarkPhoton, TestScatteringAnglePDF)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double sigma_p = pb;
	DM_Particle_Dark_Photon DM(mDM, sigma_p);
	double vDM		= 1.0e-3;
	double cosalpha = 0.2;
	obscura::Isotope target(8, 16);
	//ACT & ASSERT
	EXPECT_LT(DM.PDF_Scattering_Angle_Nucleus(-1.0, target, vDM), DM.PDF_Scattering_Angle_Nucleus(1.0, target, vDM));
	DM.Set_Low_Mass_Mode(true);
	EXPECT_DOUBLE_EQ(DM.PDF_Scattering_Angle_Nucleus(cosalpha, target, vDM), 0.5);
	EXPECT_DOUBLE_EQ(DM.PDF_Scattering_Angle_Electron(-cosalpha, vDM), 0.5);
	DM.Set_FormFactor_DM("General", keV);
	EXPECT_LT(DM.PDF_Scattering_Angle_Nucleus(-0.9, target, vDM), DM.PDF_Scattering_Angle_Nucleus(0.9, target, vDM));
	EXPECT_LT(DM.PDF_Scattering_Angle_Electron(-0.9, vDM), DM.PDF_Scattering_Angle_Electron(0.9, vDM));
}

TEST(TestDarkPhoton, TestScatteringAngleCDF)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double sigma_p = pb;
	DM_Particle_Dark_Photon DM(mDM, sigma_p);
	double vDM = 1.0e-3;
	double tol = 1e-6;

	obscura::Isotope target(8, 16);
	//ACT & ASSERT
	EXPECT_NEAR(DM.CDF_Scattering_Angle_Nucleus(-1.0, target, vDM), 0.0, tol);
	EXPECT_NEAR(DM.CDF_Scattering_Angle_Nucleus(1.0, target, vDM), 1.0, tol);
	DM.Set_Low_Mass_Mode(true);
	EXPECT_NEAR(DM.CDF_Scattering_Angle_Nucleus(-1.0, target, vDM), 0.0, tol);
	EXPECT_NEAR(DM.CDF_Scattering_Angle_Nucleus(1.0, target, vDM), 1.0, tol);
	EXPECT_NEAR(DM.CDF_Scattering_Angle_Electron(-1.0, vDM), 0.0, tol);
	EXPECT_NEAR(DM.CDF_Scattering_Angle_Electron(1.0, vDM), 1.0, tol);
}

TEST(TestDarkPhoton, TestScatteringAngleSampling)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double sigma_p = pb;
	DM_Particle_Dark_Photon DM(mDM, sigma_p);
	double vDM = 1.0e-3;
	obscura::Isotope target(8, 16);
	std::random_device rd;
	std::mt19937 PRNG(rd());

	// ACT & ASSERT
	EXPECT_LT(DM.Sample_Scattering_Angle_Electron(PRNG, vDM), 1.0);
	EXPECT_GT(DM.Sample_Scattering_Angle_Electron(PRNG, vDM), -1.0);
	EXPECT_LT(DM.Sample_Scattering_Angle_Nucleus(PRNG, target, vDM), 1.0);
	EXPECT_GT(DM.Sample_Scattering_Angle_Nucleus(PRNG, target, vDM), -1.0);
	DM.Set_Low_Mass_Mode(true);
	EXPECT_LT(DM.Sample_Scattering_Angle_Nucleus(PRNG, target, vDM), 1.0);
	EXPECT_GT(DM.Sample_Scattering_Angle_Nucleus(PRNG, target, vDM), -1.0);
	DM.Set_FormFactor_DM("General", GeV);
	EXPECT_LT(DM.Sample_Scattering_Angle_Electron(PRNG, vDM), 1.0);
	EXPECT_GT(DM.Sample_Scattering_Angle_Electron(PRNG, vDM), -1.0);
	EXPECT_LT(DM.Sample_Scattering_Angle_Nucleus(PRNG, target, vDM), 1.0);
	EXPECT_GT(DM.Sample_Scattering_Angle_Nucleus(PRNG, target, vDM), -1.0);
}

TEST(TestDarkPhoton, TestPrintSummary)
{
	// ARRANGE
	double mDM	   = 100 * MeV;
	double sigma_p = pb;
	DM_Particle_Dark_Photon DM(mDM, sigma_p);
	//ACT & ASSERT
	DM.Print_Summary();
}