#include "Scattering_Rates.hpp"

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"

namespace DaMaSCUS_SUN
{
using namespace libphysica::natural_units;
using namespace std::complex_literals;

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
	if(!use_medium_effects)
		return 1.0;
	else
	{
		double q0						 = (q * kDM * cos_theta / mDM + 0.5 * q * q / mDM);
		std::complex<double> denominator = q * q + Polarization_Tensor_L(q0, q, temperature, number_density_electron);
		return q * q * q * q / std::norm(denominator);
	}
}

double Differential_Scattering_Rate(double q, double cos_theta, double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects)
{
	double prefactor	   = electron_density * std::sqrt(2.0 / M_PI) * std::sqrt(mElectron / temperature);
	double k_1			   = DM.mass * vDM;
	double medium_function = Medium_Function(electron_density, temperature, q, DM.mass, k_1, cos_theta, use_medium_effects);
	double p1min		   = std::fabs(q / 2.0 * (1.0 + mElectron / DM.mass) + k_1 * mElectron / DM.mass * cos_theta);
	return prefactor * q * DM.dSigma_dq2_Electron(q, vDM) * vDM * vDM * medium_function * std::exp(-p1min * p1min / 2.0 / mElectron / temperature);
}

double Total_Scattering_Rate(double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects, double xi)
{

	std::function<double(double, double)> integrand = [temperature, electron_density, &DM, vDM, use_medium_effects](double q, double cos_theta) {
		return Differential_Scattering_Rate(q, cos_theta, electron_density, temperature, DM, vDM, use_medium_effects);
	};

	double vRel_max = 0.2;
	double q_min	= xi * DM.mass * vDM;
	double q_max	= 2.0 * libphysica::Reduced_Mass(DM.mass, mElectron) * vRel_max;
	double integral = libphysica::Integrate_2D(integrand, q_min, q_max, -1.0, 1.0);

	return integral;
}

double PDF_Cos_Theta(double cos_theta, double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects, double xi)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [temperature, electron_density, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate(q, cos_theta, electron_density, temperature, DM, vDM, use_medium_effects);
	};
	double vRel_max			= 0.2;
	double q_min			= xi * DM.mass * vDM;
	double q_max			= 2.0 * libphysica::Reduced_Mass(DM.mass, mElectron) * vRel_max;
	double dGamma_dcostheta = libphysica::Integrate(integrand, q_min, q_max);

	// 2. Compute total rate for normalization
	double Gamma = Total_Scattering_Rate(electron_density, temperature, DM, vDM, use_medium_effects, xi);
	return dGamma_dcostheta / Gamma;
}

double Conditional_PDF_q(double q, double cos_theta, double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects, double xi)
{
	// 1. Obtain dGamma/dcos_theta
	std::function<double(double)> integrand = [temperature, electron_density, &DM, vDM, use_medium_effects, cos_theta](double q) {
		return Differential_Scattering_Rate(q, cos_theta, electron_density, temperature, DM, vDM, use_medium_effects);
	};
	double vRel_max			= 0.2;
	double q_min			= xi * DM.mass * vDM;
	double q_max			= 2.0 * libphysica::Reduced_Mass(DM.mass, mElectron) * vRel_max;
	double dGamma_dcostheta = libphysica::Integrate(integrand, q_min, q_max);

	// 2. Compute total rate for normalization
	return Differential_Scattering_Rate(q, cos_theta, electron_density, temperature, DM, vDM, use_medium_effects) / dGamma_dcostheta;
}

double Sample_Cos_Theta(std::mt19937& PRNG, double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects, double xi)
{
	std::function<double(double)> pdf = [electron_density, temperature, &DM, vDM, use_medium_effects, xi](double cos) {
		return PDF_Cos_Theta(cos, electron_density, temperature, DM, vDM, use_medium_effects, xi);
	};
	double pdf_max = std::max(pdf(-1.0), pdf(1.0));
	return libphysica::Rejection_Sampling(pdf, -1.0, 1.0, pdf_max, PRNG);
}

double Sample_q(std::mt19937& PRNG, double cos_theta, double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects, double xi)
{
	std::function<double(double)> pdf = [cos_theta, electron_density, temperature, &DM, vDM, use_medium_effects, xi](double q) {
		return Conditional_PDF_q(q, cos_theta, electron_density, temperature, DM, vDM, use_medium_effects, xi);
	};
	double vRel_max = 0.2;
	double qMin		= xi * DM.mass * vDM;
	double qMax		= 2.0 * libphysica::Reduced_Mass(DM.mass, mElectron) * vRel_max;
	double pdf_max	= 1.1 * pdf(libphysica::Find_Maximum(pdf, qMin, qMax));
	return libphysica::Rejection_Sampling(pdf, qMin, qMax, pdf_max, PRNG);
}

}	// namespace DaMaSCUS_SUN