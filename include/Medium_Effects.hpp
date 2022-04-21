#ifndef __Medium_Effects_hpp_
#define __Medium_Effects_hpp_

#include <complex>

namespace DaMaSCUS_SUN
{

extern double Dawson_Integral(double x);
extern double Erfi(double x);

extern std::complex<double> Plasma_Dispersion_Function(double x);

}	// namespace DaMaSCUS_SUN
#endif