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
	return std::sqrt(M_PI) * std::exp(-x * x) * (1i - Erfi(x));
}

std::complex<double> Polarization_Tensor_L(double q0, double q, double temperature, double electron_number_density)
{
	double sigma_MB				= std::sqrt(temperature / mElectron);
	double plasma_frequency_sqr = Elementary_Charge * Elementary_Charge * electron_number_density / mElectron;
	double xi					= q0 / std::sqrt(2.0) / sigma_MB / q;
	double delta				= q / 2.0 / std::sqrt(2.0) / mElectron / sigma_MB;
	return plasma_frequency_sqr * mElectron / q0 * xi * (Plasma_Dispersion_Function(xi - delta) - Plasma_Dispersion_Function(xi + delta));
}

double Total_Scattering_Rate(double electron_density, double temperature, double m_DM, double sigma_e, double v_DM, double xi)
{
	double mA		 = keV;
	double qRef		 = aEM * mElectron;
	double mu_e		 = libphysica::Reduced_Mass(m_DM, mElectron);
	double k_1		 = m_DM * v_DM;
	double prefactor = electron_density * sigma_e / 2.0 / std::sqrt(2.0 * M_PI) / mu_e / mu_e * std::sqrt(mElectron / temperature);

	std::function<double(double, double)> integrand = [temperature, electron_density, k_1, m_DM, qRef, mA](double q, double cos_theta) {
		double F_DM						 = std::pow((qRef * qRef + mA * mA) / (q * q + mA * mA), 2.0);
		double q0						 = (q * k_1 * cos_theta / m_DM + 0.5 * q * q / m_DM);
		std::complex<double> denominator = q * q + Polarization_Tensor_L(q0, q, temperature, electron_density);
		double medium_factor			 = q * q * q * q / std::norm(denominator);
		double p1min					 = std::fabs(q / 2.0 * (1.0 + mElectron / m_DM) + k_1 * mElectron / m_DM * cos_theta);
		return q * F_DM * F_DM * medium_factor * std::exp(-p1min * p1min / 2.0 / mElectron / temperature);
	};

	double q_min	= xi * k_1;
	double q_max	= 2.0 * mu_e;
	double integral = libphysica::Integrate_2D(integrand, q_min, q_max, -1.0, 1.0);

	return prefactor * integral;
}

}	// namespace DaMaSCUS_SUN