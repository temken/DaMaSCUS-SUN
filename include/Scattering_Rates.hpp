#ifndef __Scattering_Rates_hpp_
#define __Scattering_Rates_hpp_

#include <complex>
#include <string>

#include "obscura/DM_Particle.hpp"

#include "Solar_Model.hpp"

namespace DaMaSCUS_SUN
{

// 1. Differential scattering rate dGamma / dq / dcos(theta)
extern double Differential_Scattering_Rate_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects);
extern double Differential_Scattering_Rate_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects);

// 2. Total scattering rate
extern double Total_Scattering_Rate_Electron(obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects, double qMin, double qMax);
extern double Total_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects, double qMin, double qMax);

// 3. Thermal average of relative speed between a particle of speed v_DM and a solar thermal target.
extern double Thermal_Averaged_Relative_Speed(double temperature, double mass_target, double v_DM);

// PDFs and sampling functions of cos(theta) and q (WILL BE MOVED TO ANOTHER FILE)
extern double PDF_Cos_Theta_Electron(double cos_theta, obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects, double qMin, double qMax);
extern double PDF_Cos_Theta_Nucleus(double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects, double qMin, double qMax);

extern double PDF_q_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects, double qMin, double qMax);
extern double PDF_q_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects, double qMin, double qMax);

extern double CDF_q_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects, double qMin, double qMax);
extern double CDF_q_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects, double qMin, double qMax);

extern double Sample_Cos_Theta_Electron(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects, double qMin, double qMax);
extern double Sample_Cos_Theta_Nucleus(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects, double qMin, double qMax);

extern double Sample_q_Electron(std::mt19937& PRNG, double cos_theta, obscura::DM_Particle& DM, double vDM, Plasma& plasma, bool use_medium_effects, double qMin, double qMax);
extern double Sample_q_Nucleus(std::mt19937& PRNG, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, Plasma& plasma, bool use_medium_effects, double qMin, double qMax);

}	// namespace DaMaSCUS_SUN
#endif