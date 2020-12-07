#ifndef __Dark_Photon_hpp_
#define __Dark_Photon_hpp_

// Headers from obscura
#include "DM_Particle.hpp"
#include "Target_Nucleus.hpp"

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

	// Differential cross sections with nuclear isotopes, elements, and electrons
	virtual double dSigma_dq2_Nucleus(double q, const obscura::Isotope& target, double vDM) const override;
	virtual double dSigma_dq2_Electron(double q, double vDM) const override;

	// Total cross sections with nuclear isotopes, elements, and electrons
	virtual double Sigma_Proton() const override;
	virtual double Sigma_Electron() const override;
	virtual double Sigma_Nucleus(const obscura::Isotope& target, double vDM) const override;

	virtual void Print_Summary(int MPI_rank = 0) const override;
};
}	// namespace DaMaSCUS_SUN
#endif