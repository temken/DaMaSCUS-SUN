#include "Medium_Effects.hpp"

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"

namespace DaMaSCUS_SUN
{
using namespace libphysica::natural_units;
using namespace std::complex_literals;

double Dawson_Integral(double x)
{
	double ans;
	static const double H = 0.4, A1 = 2.0 / 3.0, A2 = 0.4, A3 = 2.0 / 7.0;

	if(std::fabs(x) < 0.2)
	{
		double x2 = x * x;
		ans		  = x * (1.0 - A1 * x2 * (1.0 - A2 * x2 * (1.0 - A3 * x2)));
	}
	else
	{
		static const int NMAX = 6;
		static std::vector<double> c(NMAX);
		for(int i = 0; i < NMAX; i++)
			c[i] = exp(-(2.0 * i + 1.0) * (2.0 * i + 1.0) * H * H);

		double xx  = std::fabs(x);
		int n0	   = 2 * int(0.5 * xx / H + 0.5);
		double xp  = xx - n0 * H;
		double e1  = std::exp(2.0 * xp * H);
		double e2  = e1 * e1;
		double d1  = n0 + 1;
		double d2  = d1 - 2.0;
		double sum = 0.0;
		for(int i = 0; i < NMAX; i++, d1 += 2.0, d2 -= 2.0, e1 *= e2)
			sum += c[i] * (e1 / d1 + 1.0 / (d2 * e1));
		ans = 0.5641895835 * libphysica::Sign(std::exp(-xp * xp), x) * sum;
	}
	return ans;
}

double Erfi(double x)
{
	return 2.0 / std::sqrt(M_PI) * std::exp(x * x) * Dawson_Integral(x);
}
std::complex<double> Plasma_Dispersion_Function(double x)
{
	double expmx2_erfi = 2.0 / std::sqrt(M_PI) * Dawson_Integral(x);   // remove exp(-x^2) exp(x^2) to avoid inf value
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
	double k_1 = DM.mass * vDM;

	double medium_function = Medium_Function(electron_density, temperature, q, DM.mass, k_1, cos_theta, use_medium_effects);
	double p1min		   = std::fabs(q / 2.0 * (1.0 + mElectron / DM.mass) + k_1 * mElectron / DM.mass * cos_theta);
	return q * DM.dSigma_dq2_Electron(q, vDM) * vDM * vDM * medium_function * std::exp(-p1min * p1min / 2.0 / mElectron / temperature);
}

double Total_Scattering_Rate(double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects, double xi)
{
	double prefactor = electron_density * std::sqrt(2.0 / M_PI) * std::sqrt(mElectron / temperature);

	std::function<double(double, double)> integrand = [temperature, electron_density, &DM, vDM, use_medium_effects](double q, double cos_theta) {
		return Differential_Scattering_Rate(q, cos_theta, electron_density, temperature, DM, vDM, use_medium_effects);
	};

	double q_min	= xi * DM.mass * vDM;
	double q_max	= 2.0 * libphysica::Reduced_Mass(DM.mass, mElectron);
	double integral = libphysica::Integrate_2D(integrand, q_min, q_max, -1.0, 1.0);

	return prefactor * integral;
}

}	// namespace DaMaSCUS_SUN