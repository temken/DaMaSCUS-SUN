#include "Dark_Photon.hpp"

#include <iostream>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

DM_Particle_Dark_Photon::DM_Particle_Dark_Photon()
: DM_Particle(), alpha_dark(aEM), q_reference(aEM * mElectron), FF_DM("Contact"), m_dark_photon(GeV)
{
	using_cross_section = true;
	DD_use_eta_function = true;
	Set_Sigma_Proton(1.0e-40 * cm * cm);
}

DM_Particle_Dark_Photon::DM_Particle_Dark_Photon(double mDM)
: DM_Particle(mDM), alpha_dark(aEM), q_reference(aEM * mElectron), FF_DM("Contact"), m_dark_photon(GeV)
{
	using_cross_section = true;
	DD_use_eta_function = true;
	Set_Sigma_Proton(1.0e-40 * cm * cm);
}

DM_Particle_Dark_Photon::DM_Particle_Dark_Photon(double mDM, double sigma_p)
: DM_Particle(mDM), alpha_dark(aEM), q_reference(aEM * mElectron), FF_DM("Contact"), m_dark_photon(GeV)
{
	using_cross_section = true;
	DD_use_eta_function = true;
	Set_Sigma_Proton(sigma_p);
}

double DM_Particle_Dark_Photon::FormFactor2_DM(double q) const
{
	double FF;
	if(FF_DM == "Contact")
		FF = 1.0;
	else if(FF_DM == "General")
		FF = (q_reference * q_reference + m_dark_photon * m_dark_photon) / (q * q + m_dark_photon * m_dark_photon);
	else if(FF_DM == "Long-Range")
		FF = q_reference * q_reference / q / q;
	else if(FF_DM == "Electric-Dipole")
		FF = q_reference / q;
	else
	{
		std::cerr << "Error in obscura::DM_Particle_Dark_Photon::FormFactor2_DM(): Form factor " << FF_DM << "not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	return FF * FF;
}

void DM_Particle_Dark_Photon::Set_Mass(double mDM)
{
	double sigma_e = Sigma_Electron();
	mass		   = mDM;
	Set_Sigma_Electron(sigma_e);
}

void DM_Particle_Dark_Photon::Set_Epsilon(double e)
{
	epsilon = e;
}

double DM_Particle_Dark_Photon::Get_Epsilon()
{
	return epsilon;
}

// Dark matter form factor
void DM_Particle_Dark_Photon::Set_FormFactor_DM(std::string ff, double mMed)
{
	if(ff == "Contact" || ff == "Electric-Dipole" || ff == "Long-Range" || ff == "General")
		FF_DM = ff;
	else
	{
		std::cerr << "Error in obscura::DM_Particle_Dark_Photon::Set_FormFactor_DM(): Form factor " << ff << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(FF_DM == "General" && mMed > 0.0)
		m_dark_photon = mMed;
	else if(FF_DM == "Long-Range")
		m_dark_photon = 0.0;
}

void DM_Particle_Dark_Photon::Set_Dark_Photon_Mass(double m)
{
	if(FF_DM != "Long-Range")
		m_dark_photon = m;
}

// Primary interaction parameter, such as a coupling constant or cross section
double DM_Particle_Dark_Photon::Get_Interaction_Parameter(std::string target) const
{
	if(target == "Nuclei")
		return Sigma_Proton();
	else if(target == "Electrons")
		return Sigma_Electron();
	else
	{
		std::cerr << "Error in obscura::DM_Particle_Dark_Photon::Get_Interaction_Parameter(std::string): Target " << target << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void DM_Particle_Dark_Photon::Set_Interaction_Parameter(double par, std::string target)
{
	if(target == "Nuclei")
		Set_Sigma_Proton(par);
	else if(target == "Electrons")
		Set_Sigma_Electron(par);
	else
	{
		std::cerr << "Error in obscura::DM_Particle_Dark_Photon::Get_Interaction_Parameter(std::string): Target " << target << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void DM_Particle_Dark_Photon::Set_Sigma_Proton(double sigma)
{
	epsilon = (q_reference * q_reference + m_dark_photon * m_dark_photon) / 4.0 / libphysica::Reduced_Mass(mass, mProton) * sqrt(sigma / M_PI / aEM / alpha_dark);
}

void DM_Particle_Dark_Photon::Set_Sigma_Electron(double sigma)
{
	epsilon = (q_reference * q_reference + m_dark_photon * m_dark_photon) / 4.0 / libphysica::Reduced_Mass(mass, mElectron) * sqrt(sigma / M_PI / aEM / alpha_dark);
}

// Differential cross sections for nuclear targets
double DM_Particle_Dark_Photon::dSigma_dq2_Nucleus(double q, const obscura::Isotope& target, double vDM, double r) const
{
	double nuclear_form_factor = (low_mass) ? 1.0 : target.Helm_Form_Factor(q);
	double mu				   = libphysica::Reduced_Mass(mass, mProton);
	return Sigma_Proton() / 4.0 / mu / mu / vDM / vDM * FormFactor2_DM(q) * nuclear_form_factor * nuclear_form_factor * target.Z * target.Z;
}

// Differential cross section for electron targets
double DM_Particle_Dark_Photon::dSigma_dq2_Electron(double q, double vDM, double r) const
{
	double mu = libphysica::Reduced_Mass(mass, mElectron);
	return Sigma_Electron() / 4.0 / mu / mu / vDM / vDM * FormFactor2_DM(q);
}

double DM_Particle_Dark_Photon::d2Sigma_dq2_dEe_Ionization(double q, double Ee, double vDM, obscura::Atomic_Electron& shell) const
{
	return 1.0 / 4.0 / Ee * dSigma_dq2_Electron(q, vDM) * shell.Ionization_Form_Factor(q, Ee);
}

double DM_Particle_Dark_Photon::d2Sigma_dq2_dEe_Crystal(double q, double Ee, double vDM, obscura::Crystal& crystal) const
{
	return 2.0 * aEM * mElectron * mElectron / q / q / q * dSigma_dq2_Electron(q, vDM) * crystal.Crystal_Form_Factor(q, Ee);
}

// Total cross sections with nuclear isotopes, elements, and electrons
double DM_Particle_Dark_Photon::Sigma_Proton() const
{
	double mu = libphysica::Reduced_Mass(mass, mProton);
	return 16.0 * M_PI * aEM * alpha_dark * epsilon * epsilon * mu * mu / pow((q_reference * q_reference + m_dark_photon * m_dark_photon), 2.0);
}

double DM_Particle_Dark_Photon::Sigma_Electron() const
{
	double mu = libphysica::Reduced_Mass(mass, mElectron);
	return 16.0 * M_PI * aEM * alpha_dark * epsilon * epsilon * mu * mu / pow((q_reference * q_reference + m_dark_photon * m_dark_photon), 2.0);
}

double DM_Particle_Dark_Photon::Sigma_Total_Nucleus(const obscura::Isotope& target, double vDM, double r)
{
	double sigmatot = 0.0;
	if(FF_DM != "Contact" && FF_DM != "General")
	{
		std::cerr << "Error in obscura::DM_Particle_Dark_Photon::Sigma_Nucleus(): Divergence in the IR." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(!low_mass)
		sigmatot = Sigma_Total_Nucleus_Base(target, vDM, r);
	else
	{
		double mu_p = libphysica::Reduced_Mass(mass, mProton);
		double mu_N = libphysica::Reduced_Mass(mass, target.mass);
		sigmatot	= Sigma_Proton() * mu_N * mu_N / mu_p / mu_p * target.Z * target.Z;
		if(FF_DM == "General")
		{
			double q2max = 4.0 * pow(mu_N * vDM, 2.0);
			sigmatot *= pow(q_reference * q_reference + m_dark_photon * m_dark_photon, 2.0) / m_dark_photon / m_dark_photon / (m_dark_photon * m_dark_photon + q2max);
		}
	}
	return sigmatot;
}

double DM_Particle_Dark_Photon::Sigma_Total_Electron(double vDM, double r)
{
	double sigmatot = 0.0;
	if(FF_DM != "Contact" && FF_DM != "General")
	{
		std::cerr << "Error in DM_Particle_Dark_Photon::Sigma_Total_Electron(): Divergence in the IR." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		sigmatot = Sigma_Electron();
		if(FF_DM == "General")
		{
			double q2max = 4.0 * pow(libphysica::Reduced_Mass(mass, mElectron) * vDM, 2.0);
			sigmatot *= pow(q_reference * q_reference + m_dark_photon * m_dark_photon, 2.0) / m_dark_photon / m_dark_photon / (m_dark_photon * m_dark_photon + q2max);
		}
	}
	return sigmatot;
}

// Scattering angle functions
double DM_Particle_Dark_Photon::PDF_Scattering_Angle_Nucleus(double cos_alpha, const obscura::Isotope& target, double vDM, double r)
{
	if(FF_DM != "Contact" && FF_DM != "General")
	{
		std::cerr << "Error in DM_Particle_Dark_Photon::PDF_Scattering_Angle_Nucleus(): Divergence in the IR." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(!low_mass)
		return PDF_Scattering_Angle_Nucleus_Base(cos_alpha, target, vDM);
	else if(FF_DM == "Contact")
		return 0.5;
	else
	{
		double m2	 = m_dark_photon * m_dark_photon;
		double q2max = 4.0 * pow(libphysica::Reduced_Mass(mass, target.mass) * vDM, 2.0);
		return 2.0 * m2 * (m2 + q2max) / pow(2 * m2 + q2max * (1.0 - cos_alpha), 2.0);
	}
}
double DM_Particle_Dark_Photon::PDF_Scattering_Angle_Electron(double cos_alpha, double vDM, double r)
{
	if(FF_DM != "Contact" && FF_DM != "General")
	{
		std::cerr << "Error in DM_Particle_Dark_Photon::PDF_Scattering_Angle_Electron(): Divergence in the IR." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(FF_DM == "Contact")
		return 0.5;
	else
	{
		double m2	 = m_dark_photon * m_dark_photon;
		double q2max = 4.0 * pow(libphysica::Reduced_Mass(mass, mElectron) * vDM, 2.0);
		return 2.0 * m2 * (m2 + q2max) / pow(2 * m2 + q2max * (1.0 - cos_alpha), 2.0);
	}
}
double DM_Particle_Dark_Photon::CDF_Scattering_Angle_Nucleus(double cos_alpha, const obscura::Isotope& target, double vDM, double r)
{
	if(FF_DM != "Contact" && FF_DM != "General")
	{
		std::cerr << "Error in DM_Particle_Dark_Photon::CDF_Scattering_Angle_Nucleus(): Divergence in the IR." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(!low_mass)
		return CDF_Scattering_Angle_Nucleus_Base(cos_alpha, target, vDM);
	else if(FF_DM == "Contact")
		return (1.0 + cos_alpha) / 2.0;
	else
	{
		double m2	 = m_dark_photon * m_dark_photon;
		double q2max = 4.0 * pow(libphysica::Reduced_Mass(mass, target.mass) * vDM, 2.0);
		return (1.0 + cos_alpha) * m2 / (2.0 * m2 + q2max * (1.0 - cos_alpha));
	}
}

double DM_Particle_Dark_Photon::CDF_Scattering_Angle_Electron(double cos_alpha, double vDM, double r)
{
	if(FF_DM != "Contact" && FF_DM != "General")
	{
		std::cerr << "Error in DM_Particle_Dark_Photon::CDF_Scattering_Angle_Electron(): Divergence in the IR." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(FF_DM == "Contact")
		return (1.0 + cos_alpha) / 2.0;
	else
	{
		double m2	 = m_dark_photon * m_dark_photon;
		double q2max = 4.0 * pow(libphysica::Reduced_Mass(mass, mElectron) * vDM, 2.0);
		return (1.0 + cos_alpha) * m2 / (2.0 * m2 + q2max * (1.0 - cos_alpha));
	}
}

double DM_Particle_Dark_Photon::Sample_Scattering_Angle_Nucleus(std::mt19937& PRNG, const obscura::Isotope& target, double vDM, double r)
{
	if(FF_DM != "Contact" && FF_DM != "General")
	{
		std::cerr << "Error in DM_Particle_Dark_Photon::Sample_Scattering_Angle_Nucleus(): Divergence in the IR." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(!low_mass)
		return Sample_Scattering_Angle_Nucleus_Base(PRNG, target, vDM, r);
	else if(FF_DM == "Contact")
	{
		double xi = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
		return 2.0 * xi - 1.0;
	}
	else
	{
		double xi	 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
		double m2	 = m_dark_photon * m_dark_photon;
		double q2max = 4.0 * pow(libphysica::Reduced_Mass(mass, target.mass) * vDM, 2.0);
		return (m2 * (2.0 * xi - 1.0) + q2max * xi) / (m2 + q2max * xi);
	}
}

double DM_Particle_Dark_Photon::Sample_Scattering_Angle_Electron(std::mt19937& PRNG, double vDM, double r)
{
	double xi = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
	if(FF_DM != "Contact" && FF_DM != "General")
	{
		std::cerr << "Error in DM_Particle_Dark_Photon::Sample_Scattering_Angle_Electron(): Divergence in the IR." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(FF_DM == "Contact")
		return 2.0 * xi - 1.0;
	else
	{
		double m2	 = m_dark_photon * m_dark_photon;
		double q2max = 4.0 * pow(libphysica::Reduced_Mass(mass, mElectron) * vDM, 2.0);
		return (m2 * (2.0 * xi - 1.0) + q2max * xi) / (m2 + q2max * xi);
	}
}

void DM_Particle_Dark_Photon::Print_Summary(int MPI_rank) const
{
	if(MPI_rank == 0)
	{
		Print_Summary_Base();
		double massunit			= (m_dark_photon < keV) ? eV : ((m_dark_photon < MeV) ? keV : ((m_dark_photon < GeV) ? MeV : GeV));
		std::string massunitstr = (m_dark_photon < keV) ? "[eV]" : ((m_dark_photon < MeV) ? "[keV]" : ((m_dark_photon < GeV) ? "[MeV]" : "[GeV]"));
		std::cout << std::endl
				  << "\tInteraction:\t\tDark photon (DP)" << std::endl
				  << std::endl;
		std::cout << "\tInteraction type:\t" << FF_DM << std::endl
				  << "\tKinetic mixing:\t\t" << libphysica::Round(epsilon) << std::endl
				  << "\tGauge coupling a_D:\t" << libphysica::Round(alpha_dark) << std::endl
				  << "\tDark photon mass" << massunitstr << ":\t" << libphysica::Round(In_Units(m_dark_photon, massunit)) << std::endl
				  << "\tSigma_P[cm^2]:\t\t" << libphysica::Round(In_Units(Sigma_Proton(), cm * cm)) << std::endl
				  << "\tSigma_N[cm^2]:\t\t" << libphysica::Round(In_Units(Sigma_Neutron(), cm * cm)) << std::endl
				  << "\tSigma_E[cm^2]:\t\t" << libphysica::Round(In_Units(Sigma_Electron(), cm * cm)) << std::endl
				  << std::endl;

		std::cout
			<< "----------------------------------------" << std::endl;
	}
}
}	// namespace DaMaSCUS_SUN
