#ifndef __Medium_Effects_hpp_
#define __Medium_Effects_hpp_

#include <complex>

namespace DaMaSCUS_SUN
{

extern double Dawson_Integral(double x);
extern double Erfi(double x);

extern std::complex<double> Plasma_Dispersion_Function(double x);

extern std::complex<double> Polarization_Tensor_L(double q0, double q, double temperature, double electron_number_density);

}	// namespace DaMaSCUS_SUN
#endif