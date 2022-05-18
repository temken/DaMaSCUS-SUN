#include "Scattering_Rates.hpp"

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

namespace DaMaSCUS_SUN
{
using namespace libphysica::natural_units;

// 1. Differential scattering rate dGamma / dq / dcos(theta)
double Differential_Scattering_Rate_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects)
{
	double prefactor = electron_density * std::sqrt(2.0 / M_PI) * std::sqrt(mElectron / temperature);
	double k_1		 = DM.mass * vDM;
	// To DO: Medium function only contains nuclear density, which is wrong here!! Needs to be fixed.
	// double medium_function = use_medium_effects ? Medium_Function(electron_density, temperature, q, DM.mass, k_1, cos_theta, use_medium_effects) : 1.0;
	double medium_function = 1.0;
	double p1min		   = std::fabs(q / 2.0 * (1.0 + mElectron / DM.mass) + k_1 * mElectron / DM.mass * cos_theta);
	if(vDM == 0.0)	 // cancels in next expression
		vDM = 1.0;
	return prefactor * q * DM.dSigma_dq2_Electron(q, vDM) * vDM * vDM * medium_function * std::exp(-p1min * p1min / 2.0 / mElectron / temperature);
}

double Differential_Scattering_Rate_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects)
{
	double mNucleus	 = target.mass;
	double prefactor = nucleus_density * std::sqrt(2.0 / M_PI) * std::sqrt(mNucleus / temperature);
	double k_1		 = DM.mass * vDM;
	// To DO: Medium function only contains nuclear density, which is wrong here!! Needs to be fixed.
	// double medium_function	= use_medium_effects ? Medium_Function(electron_density, temperature, q, DM.mass, k_1, cos_theta, use_medium_effects) : 1.0;
	double medium_function = 1.0;
	double p1min		   = std::fabs(q / 2.0 * (1.0 + mNucleus / DM.mass) + k_1 * mNucleus / DM.mass * cos_theta);
	if(vDM == 0.0)	 // cancels in next expression
		vDM = 1.0;
	return prefactor * q * DM.dSigma_dq2_Nucleus(q, target, vDM) * vDM * vDM * medium_function * std::exp(-p1min * p1min / 2.0 / mNucleus / temperature);
}

// 2. Total scattering rate
double Total_Scattering_Rate_Electron(obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects, double qMin, double qMax)
{
	if(DM.Is_Sigma_Total_V_Dependent() || use_medium_effects || qMin > 0.0)
	{
		std::function<double(double, double)> integrand = [temperature, electron_density, &DM, vDM, use_medium_effects](double q, double cos_theta) {
			return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, electron_density, temperature, use_medium_effects);
		};

		double integral = libphysica::Integrate_2D(integrand, qMin, qMax, -1.0, 1.0);

		return integral;
	}
	else
	{
		double v_rel = Thermal_Averaged_Relative_Speed(temperature, mElectron, vDM);
		return electron_density * DM.Sigma_Total_Electron(vDM) * v_rel;
	}
}

double Total_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects, double qMin, double qMax)
{
	double mNucleus = target.mass;
	if(DM.Is_Sigma_Total_V_Dependent() || use_medium_effects || qMin > 0.0)
	{
		std::function<double(double, double)> integrand = [&target, temperature, nucleus_density, &DM, vDM, use_medium_effects](double q, double cos_theta) {
			return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, target, nucleus_density, temperature, use_medium_effects);
		};

		double integral = libphysica::Integrate_2D(integrand, qMin, qMax, -1.0, 1.0);

		return integral;
	}
	else
	{
		double v_rel = Thermal_Averaged_Relative_Speed(temperature, mNucleus, vDM);
		return nucleus_density * DM.Sigma_Total_Nucleus(target, vDM) * v_rel;
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
double PDF_Cos_Theta_Electron(double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects, double qMin, double qMax)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [temperature, electron_density, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, electron_density, temperature, use_medium_effects);
	};
	double dGamma_dcostheta = libphysica::Integrate(integrand, qMin, qMax);

	// 2. Compute total rate for normalization
	double Gamma = Total_Scattering_Rate_Electron(DM, vDM, electron_density, temperature, use_medium_effects, qMin, qMax);
	return dGamma_dcostheta / Gamma;
}

double PDF_Cos_Theta_Nucleus(double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects, double qMin, double qMax)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [&target, temperature, nucleus_density, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, target, nucleus_density, temperature, use_medium_effects);
	};
	double dGamma_dcostheta = libphysica::Integrate(integrand, qMin, qMax);

	// 2. Compute total rate for normalization
	double Gamma = Total_Scattering_Rate_Nucleus(DM, vDM, target, nucleus_density, temperature, use_medium_effects, qMin, qMax);
	return dGamma_dcostheta / Gamma;
}

double PDF_q_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects, double qMin, double qMax)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [temperature, electron_density, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, electron_density, temperature, use_medium_effects);
	};
	double dGamma_dcostheta = libphysica::Integrate(integrand, qMin, qMax);

	// 2. Compute total rate for normalization
	return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, electron_density, temperature, use_medium_effects) / dGamma_dcostheta;
}

double PDF_q_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects, double qMin, double qMax)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [&target, temperature, nucleus_density, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, target, nucleus_density, temperature, use_medium_effects);
	};
	double dGamma_dcostheta = libphysica::Integrate(integrand, qMin, qMax);

	// 2. Compute total rate for normalization
	return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, target, nucleus_density, temperature, use_medium_effects) / dGamma_dcostheta;
}

double CDF_q_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects, double qMin, double qMax)
{
	std::function<double(double, double)> integrand = [temperature, electron_density, &DM, vDM, use_medium_effects](double q, double cos_theta) {
		return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, electron_density, temperature, use_medium_effects);
	};

	double integral	  = libphysica::Integrate_2D(integrand, qMin, q, -1.0, 1.0);
	double integral_2 = libphysica::Integrate_2D(integrand, qMin, qMax, -1.0, 1.0);

	return integral / integral_2;
}

double CDF_q_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects, double qMin, double qMax)
{
	std::function<double(double, double)> integrand = [temperature, &target, nucleus_density, &DM, vDM, use_medium_effects](double q, double cos_theta) {
		return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, target, nucleus_density, temperature, use_medium_effects);
	};

	double integral	  = libphysica::Integrate_2D(integrand, qMin, q, -1.0, 1.0);
	double integral_2 = libphysica::Integrate_2D(integrand, qMin, qMax, -1.0, 1.0);

	return integral / integral_2;
}

double Sample_Cos_Theta_Electron(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects, double qMin, double qMax)
{
	std::function<double(double)> pdf = [electron_density, temperature, &DM, vDM, use_medium_effects, qMin, qMax](double cos) {
		return PDF_Cos_Theta_Electron(cos, DM, vDM, electron_density, temperature, use_medium_effects, qMin, qMax);
	};
	double pdf_max = std::max(pdf(-1.0), pdf(1.0));
	return libphysica::Rejection_Sampling(pdf, -1.0, 1.0, pdf_max, PRNG);
}

double Sample_Cos_Theta_Nucleus(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects, double qMin, double qMax)
{
	std::function<double(double)> pdf = [&target, nucleus_density, temperature, &DM, vDM, use_medium_effects, qMin, qMax](double cos) {
		return PDF_Cos_Theta_Nucleus(cos, DM, vDM, target, nucleus_density, temperature, use_medium_effects, qMin, qMax);
	};
	double pdf_max = std::max(pdf(-1.0), pdf(1.0));
	return libphysica::Rejection_Sampling(pdf, -1.0, 1.0, pdf_max, PRNG);
}

double Sample_q_Electron(std::mt19937& PRNG, double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects, double qMin, double qMax)
{
	std::function<double(double)> cdf = [cos_theta, electron_density, temperature, &DM, vDM, use_medium_effects, qMin, qMax](double q) {
		return CDF_q_Electron(q, cos_theta, DM, vDM, electron_density, temperature, use_medium_effects, qMin, qMax);
	};
	return libphysica::Inverse_Transform_Sampling(cdf, qMin, qMax, PRNG);
}

double Sample_q_Nucleus(std::mt19937& PRNG, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects, double qMin, double qMax)
{
	std::function<double(double)> cdf = [&target, cos_theta, nucleus_density, temperature, &DM, vDM, use_medium_effects, qMin, qMax](double q) {
		return CDF_q_Nucleus(q, cos_theta, DM, vDM, target, nucleus_density, temperature, use_medium_effects, qMin, qMax);
	};
	return libphysica::Inverse_Transform_Sampling(cdf, qMin, qMax, PRNG);
}

}	// namespace DaMaSCUS_SUN