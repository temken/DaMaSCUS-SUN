#ifndef __Scattering_Rates_hpp_
#define __Scattering_Rates_hpp_

#include <complex>
#include <string>

#include "obscura/DM_Particle.hpp"

namespace DaMaSCUS_SUN
{

// 1. Differential scattering rate dGamma / dq / dcos(theta)
extern double Differential_Scattering_Rate_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, double target_density, double temperature, bool use_medium_effects = false);
extern double Differential_Scattering_Rate_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double target_density, double temperature, bool use_medium_effects = false);

// 2. Total scattering rate
extern double Total_Scattering_Rate_Electron(obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects = false, double zeta = 0.0);
extern double Total_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects = false, double zeta = 0.0);

// 3. Thermal average of relative speed between a particle of speed v_DM and a solar thermal target.
extern double Thermal_Averaged_Relative_Speed(double temperature, double mass_target, double v_DM);

// 4. Medium effects
extern std::complex<double> Plasma_Dispersion_Function(double x);
extern std::complex<double> Polarization_Tensor_L(double q0, double q, double temperature, double electron_number_density);
extern double Medium_Function(double number_density_electron, double temperature, double q, double mDM, double kDM, double cos_theta, bool use_medium_effects = false);

// PDFs and sampling functions of cos(theta) and q (WILL BE MOVED TO ANOTHER FILE)
extern double PDF_Cos_Theta_Electron(double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects = false, double zeta = 0.0);
extern double PDF_Cos_Theta_Nucleus(double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects = false, double zeta = 0.0);

extern double PDF_q_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects = false, double zeta = 0.0);
extern double PDF_q_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects = false, double zeta = 0.0);

extern double CDF_q_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects = false, double zeta = 0.0);
extern double CDF_q_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects = false, double zeta = 0.0);

extern double Sample_Cos_Theta_Electron(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects = false, double zeta = 0.0);
extern double Sample_Cos_Theta_Nucleus(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects = false, double zeta = 0.0);

extern double Sample_q_Electron(std::mt19937& PRNG, double cos_theta, obscura::DM_Particle& DM, double vDM, double electron_density, double temperature, bool use_medium_effects = false, double zeta = 0.0);
extern double Sample_q_Nucleus(std::mt19937& PRNG, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& target, double nucleus_density, double temperature, bool use_medium_effects = false, double zeta = 0.0);

}	// namespace DaMaSCUS_SUN
#endif