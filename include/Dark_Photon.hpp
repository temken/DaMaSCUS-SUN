#ifndef __Dark_Photon_hpp_
#define __Dark_Photon_hpp_

#include "obscura/DM_Particle.hpp"
#include "obscura/Target_Nucleus.hpp"

#include "Solar_Model.hpp"

namespace DaMaSCUS_SUN
{

// Dark photon model with kinetic mixing
class DM_Particle_Dark_Photon : public obscura::DM_Particle
{
  private:
	double epsilon, alpha_dark;

	// Dark matter form factor
	double q_reference;
	std::string FF_DM;
	double m_dark_photon;
	double FormFactor2_DM(double q) const;

	Solar_Model SSM;

  public:
	// Constructors
	DM_Particle_Dark_Photon();
	explicit DM_Particle_Dark_Photon(double mDM);
	explicit DM_Particle_Dark_Photon(double mDM, double sigma_p);

	virtual void Set_Mass(double mDM) override;

	void Set_Epsilon(double e);
	double Get_Epsilon();

	// Dark matter form factor
	void Set_FormFactor_DM(std::string ff, double mMed = -1.0);
	void Set_Dark_Photon_Mass(double m);

	// Primary interaction parameter, such as a coupling constant or cross section
	virtual double Get_Interaction_Parameter(std::string target) const override;
	virtual void Set_Interaction_Parameter(double par, std::string target) override;

	virtual void Set_Sigma_Proton(double sigma) override;
	virtual void Set_Sigma_Electron(double sigma) override;

	// Differential cross sections for nuclear targets
	virtual double dSigma_dq2_Nucleus(double q, const obscura::Isotope& target, double vDM, double r = -1.0) const override;

	// Differential cross section for electron targets
	virtual double dSigma_dq2_Electron(double q, double vDM, double r = -1.0) const override;
	virtual double d2Sigma_dq2_dEe_Ionization(double q, double Ee, double vDM, obscura::Atomic_Electron& shell) const override;
	virtual double d2Sigma_dq2_dEe_Crystal(double q, double Ee, double vDM, obscura::Crystal& crystal) const override;

	// Total cross sections with nuclear isotopes, elements, and electrons
	virtual double Sigma_Proton() const override;
	virtual double Sigma_Electron() const override;

	virtual bool Is_Sigma_Total_V_Dependent() const override;
	virtual double Sigma_Total_Nucleus(const obscura::Isotope& target, double vDM, double r = -1.0) override;
	virtual double Sigma_Total_Electron(double vDM, double r = -1.0) override;

	// Scattering angle functions
	virtual double PDF_Scattering_Angle_Nucleus(double cos_alpha, const obscura::Isotope& target, double vDM, double r = -1.0) override;
	virtual double PDF_Scattering_Angle_Electron(double cos_alpha, double vDM, double r = -1.0) override;
	virtual double CDF_Scattering_Angle_Nucleus(double cos_alpha, const obscura::Isotope& target, double vDM, double r = -1.0) override;
	virtual double CDF_Scattering_Angle_Electron(double cos_alpha, double vDM, double r = -1.0) override;
	virtual double Sample_Scattering_Angle_Nucleus(std::mt19937& PRNG, const obscura::Isotope& target, double vDM, double r = -1.0) override;
	virtual double Sample_Scattering_Angle_Electron(std::mt19937& PRNG, double vDM, double r = -1.0) override;

	virtual void Print_Summary(int MPI_rank = 0) const override;
};
}	// namespace DaMaSCUS_SUN
#endif