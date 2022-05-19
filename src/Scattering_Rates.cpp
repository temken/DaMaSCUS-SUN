#include "Scattering_Rates.hpp"

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

namespace DaMaSCUS_SUN
{
using namespace libphysica::natural_units;

// 1. Differential scattering rate dGamma / dq / dcos(theta)
double Differential_Scattering_Rate_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects)
{
	double prefactor = plasma.number_density_electrons * std::sqrt(2.0 / M_PI) * std::sqrt(mElectron / plasma.temperature);
	double k_1		 = DM.mass * vDM;
	if(use_medium_effects)
	{
		double q0 = 1.0 / 2.0 / DM.mass * (q * q + 2.0 * q * k_1 * cos_theta);
		prefactor *= plasma.Form_Factor_Medium_Effects(q0, q);
	}
	double p1min = std::fabs(q / 2.0 * (1.0 + mElectron / DM.mass) + k_1 * mElectron / DM.mass * cos_theta);
	vDM			 = 1.0;	  // cancels in next expression (in case vDM == 0)
	return prefactor * q * DM.dSigma_dq2_Electron(q, vDM) * vDM * vDM * std::exp(-p1min * p1min / 2.0 / mElectron / plasma.temperature);
}

double Differential_Scattering_Rate_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects)
{
	double prefactor = nucleus_density * std::sqrt(2.0 / M_PI) * std::sqrt(nucleus.mass / plasma.temperature);
	double k_1		 = DM.mass * vDM;
	if(use_medium_effects)
	{
		double q0 = 1.0 / 2.0 / DM.mass * (q * q + 2.0 * q * k_1 * cos_theta);
		prefactor *= plasma.Form_Factor_Medium_Effects(q0, q);
	}
	double p1min = std::fabs(q / 2.0 * (1.0 + nucleus.mass / DM.mass) + k_1 * nucleus.mass / DM.mass * cos_theta);
	vDM			 = 1.0;	  // cancels in next expression (in case vDM == 0)
	return prefactor * q * DM.dSigma_dq2_Nucleus(q, nucleus, vDM) * vDM * vDM * std::exp(-p1min * p1min / 2.0 / nucleus.mass / plasma.temperature);
}

// 2. Total scattering rate
double Total_Scattering_Rate_Electron(obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects, double qMin, double qMax)
{
	if(DM.Is_Sigma_Total_V_Dependent() || use_medium_effects || qMin > 0.0)
	{
		std::function<double(double, double)> integrand = [&DM, vDM, &plasma, use_medium_effects](double q, double cos_theta) {
			return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, plasma, use_medium_effects);
		};

		double integral = libphysica::Integrate_2D(integrand, qMin, qMax, -1.0, 1.0);
		return integral;
	}
	else
	{
		double v_rel = Thermal_Averaged_Relative_Speed(plasma.temperature, mElectron, vDM);
		return plasma.number_density_electrons * DM.Sigma_Total_Electron(vDM) * v_rel;
	}
}

double Total_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects, double qMin, double qMax)
{
	if(DM.Is_Sigma_Total_V_Dependent() || use_medium_effects || qMin > 0.0)
	{
		std::function<double(double, double)> integrand = [&nucleus, nucleus_density, &DM, vDM, &plasma, use_medium_effects](double q, double cos_theta) {
			return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, nucleus, nucleus_density, plasma, use_medium_effects);
		};
		return libphysica::Integrate_2D(integrand, qMin, qMax, -1.0, 1.0);
	}
	else
	{
		double v_rel = Thermal_Averaged_Relative_Speed(plasma.temperature, nucleus.mass, vDM);
		return nucleus_density * DM.Sigma_Total_Nucleus(nucleus, vDM) * v_rel;
	}
}

// 3. Thermal average of relative speed between a particle of speed v_DM and a solar thermal target.
double Thermal_Averaged_Relative_Speed(double temperature, double mass_target, double v_DM)
{
	double kappa = sqrt(mass_target / 2.0 / temperature);
	if(v_DM < 1.0e-20)
		return 2.0 / sqrt(M_PI) / kappa;
	else
	{
		double relative_speed = (1.0 + 2.0 * pow(kappa * v_DM, 2.0)) * erf(kappa * v_DM) / 2.0 / kappa / kappa / v_DM + exp(-pow(kappa * v_DM, 2.0)) / sqrt(M_PI) / kappa;
		return relative_speed;
	}
}

// PDFs and sampling functions of cos(theta) and q (WILL BE MOVED TO ANOTHER FILE)
double PDF_Cos_Theta_Electron(double cos_theta, obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects, double qMin, double qMax)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [&plasma, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, plasma, use_medium_effects);
	};
	double dGamma_dcostheta = libphysica::Integrate(integrand, qMin, qMax);

	// 2. Compute total rate for normalization
	double Gamma = Total_Scattering_Rate_Electron(DM, vDM, plasma, use_medium_effects, qMin, qMax);
	return dGamma_dcostheta / Gamma;
}

double PDF_Cos_Theta_Nucleus(double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects, double qMin, double qMax)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [&nucleus, &plasma, nucleus_density, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, nucleus, nucleus_density, plasma, use_medium_effects);
	};
	double dGamma_dcostheta = libphysica::Integrate(integrand, qMin, qMax);

	// 2. Compute total rate for normalization
	double Gamma = Total_Scattering_Rate_Nucleus(DM, vDM, nucleus, nucleus_density, plasma, use_medium_effects, qMin, qMax);
	return dGamma_dcostheta / Gamma;
}

// Conditional PDF/CDF of q for a fixed value of cos_theta
double PDF_q_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects, double qMin, double qMax)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [&plasma, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, plasma, use_medium_effects);
	};
	double dGamma_dcostheta = libphysica::Integrate(integrand, qMin, qMax);

	return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, plasma, use_medium_effects) / dGamma_dcostheta;
}

double PDF_q_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects, double qMin, double qMax)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [&nucleus, &plasma, nucleus_density, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, nucleus, nucleus_density, plasma, use_medium_effects);
	};
	double dGamma_dcostheta = libphysica::Integrate(integrand, qMin, qMax);

	return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, nucleus, nucleus_density, plasma, use_medium_effects) / dGamma_dcostheta;
}

double CDF_q_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects, double qMin, double qMax)
{
	std::function<double(double)> integrand = [cos_theta, &plasma, &DM, vDM, use_medium_effects](double q) {
		return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, plasma, use_medium_effects);
	};

	double integral_1 = libphysica::Integrate(integrand, qMin, q);
	double integral_2 = libphysica::Integrate(integrand, qMin, qMax);

	return integral_1 / integral_2;
}

double CDF_q_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects, double qMin, double qMax)
{
	std::function<double(double)> integrand = [cos_theta, &plasma, &nucleus, nucleus_density, &DM, vDM, use_medium_effects](double q) {
		return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, nucleus, nucleus_density, plasma, use_medium_effects);
	};

	double integral_1 = libphysica::Integrate(integrand, qMin, q);
	double integral_2 = libphysica::Integrate(integrand, qMin, qMax);

	return integral_1 / integral_2;
}

double Sample_Cos_Theta_Electron(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects, double qMin, double qMax)
{
	std::function<double(double)> pdf = [&plasma, &DM, vDM, use_medium_effects, qMin, qMax](double cos) {
		return PDF_Cos_Theta_Electron(cos, DM, vDM, plasma, use_medium_effects, qMin, qMax);
	};
	double pdf_max = std::max(pdf(-1.0), pdf(1.0));
	if(pdf_max < 1.0)
		pdf_max = 1.0;
	return libphysica::Rejection_Sampling(pdf, -1.0, 1.0, pdf_max, PRNG);
}

double Sample_Cos_Theta_Nucleus(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects, double qMin, double qMax)
{
	std::function<double(double)> pdf = [&nucleus, nucleus_density, &plasma, &DM, vDM, use_medium_effects, qMin, qMax](double cos) {
		return PDF_Cos_Theta_Nucleus(cos, DM, vDM, nucleus, nucleus_density, plasma, use_medium_effects, qMin, qMax);
	};
	double pdf_max = std::max(pdf(-1.0), pdf(1.0));
	if(pdf_max < 1.0)
		pdf_max = 1.0;
	return libphysica::Rejection_Sampling(pdf, -1.0, 1.0, pdf_max, PRNG);
}

double Sample_q_Electron(std::mt19937& PRNG, double cos_theta, obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects, double qMin, double qMax)
{
	std::function<double(double)> cdf = [cos_theta, &plasma, &DM, vDM, use_medium_effects, qMin, qMax](double q) {
		return CDF_q_Electron(q, cos_theta, DM, vDM, plasma, use_medium_effects, qMin, qMax);
	};
	return libphysica::Inverse_Transform_Sampling(cdf, qMin, qMax, PRNG);
}

double Sample_q_Nucleus(std::mt19937& PRNG, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects, double qMin, double qMax)
{
	std::function<double(double)> cdf = [&nucleus, cos_theta, nucleus_density, &plasma, &DM, vDM, use_medium_effects, qMin, qMax](double q) {
		return CDF_q_Nucleus(q, cos_theta, DM, vDM, nucleus, nucleus_density, plasma, use_medium_effects, qMin, qMax);
	};
	return libphysica::Inverse_Transform_Sampling(cdf, qMin, qMax, PRNG);
}

}	// namespace DaMaSCUS_SUN