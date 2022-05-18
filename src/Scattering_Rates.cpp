#include "Scattering_Rates.hpp"

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"

namespace DaMaSCUS_SUN
{
using namespace libphysica::natural_units;
using namespace std::complex_literals;

// 1. Differential scattering rate dGamma / dq / dcos(theta)
double Differential_Scattering_Rate_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects)
{
	double prefactor = electron_density * std::sqrt(2.0 / M_PI) * std::sqrt(mElectron / temperature);
	double k_1		 = DM.mass * vDM;
	// To DO: Medium function only contains nuclear density, which is wrong here!! Needs to be fixed.
	double medium_function = use_medium_effects ? Medium_Function(electron_density, temperature, q, DM.mass, k_1, cos_theta, use_medium_effects) : 1.0;
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
	double electron_density = 10.0 * nucleus_density;
	double medium_function	= use_medium_effects ? Medium_Function(electron_density, temperature, q, DM.mass, k_1, cos_theta, use_medium_effects) : 1.0;
	double p1min			= std::fabs(q / 2.0 * (1.0 + mNucleus / DM.mass) + k_1 * mNucleus / DM.mass * cos_theta);
	if(vDM == 0.0)	 // cancels in next expression
		vDM = 1.0;
	return prefactor * q * DM.dSigma_dq2_Nucleus(q, target, vDM) * vDM * vDM * medium_function * std::exp(-p1min * p1min / 2.0 / mNucleus / temperature);
}

// 2. Total scattering rate
double Total_Scattering_Rate_Electron(obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects, double zeta)
{
	if(DM.Is_Sigma_Total_V_Dependent() || use_medium_effects || zeta > 0.0)
	{
		std::function<double(double, double)> integrand = [temperature, electron_density, &DM, vDM, use_medium_effects](double q, double cos_theta) {
			return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, electron_density, temperature, use_medium_effects);
		};

		double vRel_max = 0.2;
		double q_min	= zeta * DM.mass * vDM;
		double q_max	= 2.0 * libphysica::Reduced_Mass(DM.mass, mElectron) * vRel_max;
		double integral = libphysica::Integrate_2D(integrand, q_min, q_max, -1.0, 1.0);

		return integral;
	}
	else
	{
		double v_rel = Thermal_Averaged_Relative_Speed(temperature, mElectron, vDM);
		return electron_density * DM.Sigma_Total_Electron(vDM) * v_rel;
	}
}

double Total_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects, double zeta)
{
	double mNucleus = target.mass;
	if(DM.Is_Sigma_Total_V_Dependent() || use_medium_effects || zeta > 0.0)
	{
		std::function<double(double, double)> integrand = [&target, temperature, nucleus_density, &DM, vDM, use_medium_effects](double q, double cos_theta) {
			return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, target, nucleus_density, temperature, use_medium_effects);
		};

		double vRel_max = 0.2;
		double q_min	= zeta * DM.mass * vDM;
		double q_max	= 2.0 * libphysica::Reduced_Mass(DM.mass, mNucleus) * vRel_max;
		double integral = libphysica::Integrate_2D(integrand, q_min, q_max, -1.0, 1.0);

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

// 4. Medium effects
std::complex<double> Plasma_Dispersion_Function(double x)
{
	double expmx2_erfi = 2.0 / std::sqrt(M_PI) * libphysica::Dawson_Integral(x);   // remove exp(-x^2) exp(x^2) to avoid inf value
	return std::sqrt(M_PI) * (1i * std::exp(-x * x) - expmx2_erfi);
}

std::complex<double> Polarization_Tensor_L(double q0, double q, double temperature, double electron_number_density)
{
	bool use_vlaslov_approximation = false;

	double sigma_MB				= std::sqrt(temperature / mElectron);
	double plasma_frequency_sqr = Elementary_Charge * Elementary_Charge * electron_number_density / mElectron;
	double xi					= q0 / std::sqrt(2.0) / sigma_MB / q;
	double delta				= q / 2.0 / std::sqrt(2.0) / mElectron / sigma_MB;

	if(!use_vlaslov_approximation)
		return plasma_frequency_sqr * mElectron / q0 * xi * (Plasma_Dispersion_Function(xi - delta) - Plasma_Dispersion_Function(xi + delta));
	else
		return plasma_frequency_sqr / sigma_MB / sigma_MB * (1.0 + xi * Plasma_Dispersion_Function(xi));
}

double Medium_Function(double number_density_electron, double temperature, double q, double mDM, double kDM, double cos_theta, bool use_medium_effects)
{
	double q0						 = (q * kDM * cos_theta / mDM + 0.5 * q * q / mDM);
	std::complex<double> denominator = q * q + Polarization_Tensor_L(q0, q, temperature, number_density_electron);
	return q * q * q * q / std::norm(denominator);
}

// PDFs and sampling functions of cos(theta) and q (WILL BE MOVED TO ANOTHER FILE)
double PDF_Cos_Theta_Electron(double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects, double zeta)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [temperature, electron_density, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, electron_density, temperature, use_medium_effects);
	};
	double vRel_max			= 0.2;
	double q_min			= zeta * DM.mass * vDM;
	double q_max			= 2.0 * libphysica::Reduced_Mass(DM.mass, mElectron) * vRel_max;
	double dGamma_dcostheta = libphysica::Integrate(integrand, q_min, q_max);

	// 2. Compute total rate for normalization
	double Gamma = Total_Scattering_Rate_Electron(DM, vDM, electron_density, temperature, use_medium_effects, zeta);
	return dGamma_dcostheta / Gamma;
}

double PDF_Cos_Theta_Nucleus(double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects, double zeta)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [&target, temperature, nucleus_density, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, target, nucleus_density, temperature, use_medium_effects);
	};
	double vRel_max			= 0.2;
	double q_min			= zeta * DM.mass * vDM;
	double q_max			= 2.0 * libphysica::Reduced_Mass(DM.mass, target.mass) * vRel_max;
	double dGamma_dcostheta = libphysica::Integrate(integrand, q_min, q_max);

	// 2. Compute total rate for normalization
	double Gamma = Total_Scattering_Rate_Nucleus(DM, vDM, target, nucleus_density, temperature, use_medium_effects, zeta);
	return dGamma_dcostheta / Gamma;
}

double PDF_q_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects, double zeta)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [temperature, electron_density, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, electron_density, temperature, use_medium_effects);
	};
	double vRel_max			= 0.2;
	double q_min			= zeta * DM.mass * vDM;
	double q_max			= 2.0 * libphysica::Reduced_Mass(DM.mass, mElectron) * vRel_max;
	double dGamma_dcostheta = libphysica::Integrate(integrand, q_min, q_max);

	// 2. Compute total rate for normalization
	return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, electron_density, temperature, use_medium_effects) / dGamma_dcostheta;
}

double PDF_q_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects, double zeta)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [&target, temperature, nucleus_density, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, target, nucleus_density, temperature, use_medium_effects);
	};
	double vRel_max			= 0.2;
	double q_min			= zeta * DM.mass * vDM;
	double q_max			= 2.0 * libphysica::Reduced_Mass(DM.mass, target.mass) * vRel_max;
	double dGamma_dcostheta = libphysica::Integrate(integrand, q_min, q_max);

	// 2. Compute total rate for normalization
	return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, target, nucleus_density, temperature, use_medium_effects) / dGamma_dcostheta;
}

double CDF_q_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects, double zeta)
{
	std::function<double(double, double)> integrand = [temperature, electron_density, &DM, vDM, use_medium_effects](double q, double cos_theta) {
		return Differential_Scattering_Rate_Electron(q, cos_theta, DM, vDM, electron_density, temperature, use_medium_effects);
	};

	double vRel_max	  = 0.2;
	double q_min	  = zeta * DM.mass * vDM;
	double q_max	  = 2.0 * libphysica::Reduced_Mass(DM.mass, mElectron) * vRel_max;
	double integral	  = libphysica::Integrate_2D(integrand, q_min, q, -1.0, 1.0);
	double integral_2 = libphysica::Integrate_2D(integrand, q_min, q_max, -1.0, 1.0);

	return integral / integral_2;
}

double CDF_q_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects, double zeta)
{
	std::function<double(double, double)> integrand = [temperature, &target, nucleus_density, &DM, vDM, use_medium_effects](double q, double cos_theta) {
		return Differential_Scattering_Rate_Nucleus(q, cos_theta, DM, vDM, target, nucleus_density, temperature, use_medium_effects);
	};

	double vRel_max	  = 0.2;
	double q_min	  = zeta * DM.mass * vDM;
	double q_max	  = 2.0 * libphysica::Reduced_Mass(DM.mass, target.mass) * vRel_max;
	double integral	  = libphysica::Integrate_2D(integrand, q_min, q, -1.0, 1.0);
	double integral_2 = libphysica::Integrate_2D(integrand, q_min, q_max, -1.0, 1.0);

	return integral / integral_2;
}

double Sample_Cos_Theta_Electron(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects, double zeta)
{
	std::function<double(double)> pdf = [electron_density, temperature, &DM, vDM, use_medium_effects, zeta](double cos) {
		return PDF_Cos_Theta_Electron(cos, DM, vDM, electron_density, temperature, use_medium_effects, zeta);
	};
	double pdf_max = std::max(pdf(-1.0), pdf(1.0));
	return libphysica::Rejection_Sampling(pdf, -1.0, 1.0, pdf_max, PRNG);
}

double Sample_Cos_Theta_Nucleus(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects, double zeta)
{
	std::function<double(double)> pdf = [&target, nucleus_density, temperature, &DM, vDM, use_medium_effects, zeta](double cos) {
		return PDF_Cos_Theta_Nucleus(cos, DM, vDM, target, nucleus_density, temperature, use_medium_effects, zeta);
	};
	double pdf_max = std::max(pdf(-1.0), pdf(1.0));
	return libphysica::Rejection_Sampling(pdf, -1.0, 1.0, pdf_max, PRNG);
}

double Sample_q_Electron(std::mt19937& PRNG, double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects, double zeta)
{
	double vRel_max					  = 0.2;
	double qMin						  = zeta * DM.mass * vDM;
	double qMax						  = 2.0 * libphysica::Reduced_Mass(DM.mass, mElectron) * vRel_max;
	double xi						  = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
	std::function<double(double)> fct = [xi, cos_theta, electron_density, temperature, &DM, vDM, use_medium_effects, zeta](double q) {
		return xi - CDF_q_Electron(q, cos_theta, DM, vDM, electron_density, temperature, use_medium_effects, zeta);
	};
	double q = libphysica::Find_Root(fct, qMin, qMax, 1e-10 * (qMax - qMin));
	return q;
}

double Sample_q_Nucleus(std::mt19937& PRNG, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects, double zeta)
{
	double vRel_max					  = 0.2;
	double qMin						  = zeta * DM.mass * vDM;
	double qMax						  = 2.0 * libphysica::Reduced_Mass(DM.mass, target.mass) * vRel_max;
	double xi						  = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
	std::function<double(double)> fct = [xi, &target, cos_theta, nucleus_density, temperature, &DM, vDM, use_medium_effects, zeta](double q) {
		return xi - CDF_q_Nucleus(q, cos_theta, DM, vDM, target, nucleus_density, temperature, use_medium_effects, zeta);
	};
	double q = libphysica::Find_Root(fct, qMin, qMax, 1e-10 * (qMax - qMin));
	return q;
}

}	// namespace DaMaSCUS_SUN