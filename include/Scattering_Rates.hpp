#ifndef __Scattering_Rates_hpp_
#define __Scattering_Rates_hpp_

#include <complex>
#include <string>

#include "obscura/DM_Particle.hpp"

namespace DaMaSCUS_SUN
{

extern std::complex<double> Plasma_Dispersion_Function(double x);

extern std::complex<double> Polarization_Tensor_L(double q0, double q, double temperature, double electron_number_density);

extern double Medium_Function(double number_density_electron, double temperature, double q, double k_1, double cos_theta, bool use_medium_effects = false);

extern double Differential_Scattering_Rate(double q, double cos_theta, double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects = false);
extern double Total_Scattering_Rate(double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects = false, double xi = 0.0);

extern double PDF_Cos_Theta(double cos_theta, double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects = false, double xi = 0.0);
extern double Conditional_PDF_q(double q, double cos_theta, double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects = false, double xi = 0.0);

extern double Sample_Cos_Theta(std::mt19937& PRNG, double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects = false, double xi = 0.0);
extern double Sample_q(std::mt19937& PRNG, double cos_theta, double electron_density, double temperature, obscura::DM_Particle& DM, double vDM, bool use_medium_effects = false, double xi = 0.0);

}	// namespace DaMaSCUS_SUN
#endif