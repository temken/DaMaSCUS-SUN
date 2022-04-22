#ifndef __Medium_Effects_hpp_
#define __Medium_Effects_hpp_

#include <complex>

namespace DaMaSCUS_SUN
{

extern double Dawson_Integral(double x);
extern double Erfi(double x);

extern std::complex<double> Plasma_Dispersion_Function(double x);

extern std::complex<double> Polarization_Tensor_L(double q0, double q, double temperature, double electron_number_density);

extern double Total_Scattering_Rate(double electron_density, double temperature, double m_DM, double sigma_e, double v_DM, double xi = 0.0);

}	// namespace DaMaSCUS_SUN
#endif