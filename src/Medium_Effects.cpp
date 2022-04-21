#include "Medium_Effects.hpp"

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

}	// namespace DaMaSCUS_SUN