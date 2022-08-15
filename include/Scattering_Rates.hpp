#ifndef __Scattering_Rates_hpp_
#define __Scattering_Rates_hpp_

#include <complex>
#include <string>

#include "obscura/DM_Particle.hpp"

#include "Solar_Model.hpp"

namespace DaMaSCUS_SUN
{

extern std::complex<double> Plasma_Dispersion_Function(double x);
extern std::complex<double> Polarization_Tensor_Longitudinal(double q0, double q, double temperature, double number_density, double mass, double Z);

template <typename Container>
extern std::complex<double> Polarization_Tensor_L(double q0, double q, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei);
template <typename Container>
extern double Form_Factor_Medium_Effects(double q0, double q, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei);

// 1. Differential scattering rate dGamma / dq / dcos(theta)
template <typename Container>
extern double Differential_Scattering_Rate_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects);
template <typename Container>
extern double Differential_Scattering_Rate_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects);

// 2. Total scattering rate
template <typename Container>
extern double Total_Scattering_Rate_Electron(obscura::DM_Particle& DM, double vDM, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);
template <typename Container>
extern double Total_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);

// 3. Thermal average of relative speed between a particle of speed v_DM and a solar thermal target.
extern double Thermal_Averaged_Relative_Speed(double temperature, double mass_target, double v_DM);
extern double Maximum_Relative_Speed(double temperature, double mass_target, double vDM, double N = 5.0);
extern double Maximum_Momentum_Transfer(double mDM, double temperature, double mass_target, double vDM, double N = 5.0);

// PDFs and CDFs and sampling functions of cos(theta) and q
template <typename Container>
extern double PDF_Cos_Theta_Electron(double cos_theta, obscura::DM_Particle& DM, double vDM, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);
template <typename Container>
extern double PDF_Cos_Theta_Nucleus(double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);

template <typename Container>
extern double CDF_Cos_Theta_Electron(double cos_theta, obscura::DM_Particle& DM, double vDM, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);
template <typename Container>
extern double CDF_Cos_Theta_Nucleus(double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);

// Conditional PDF/CDF of q for a fixed value of cos_theta
template <typename Container>
extern double PDF_q_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);
template <typename Container>
extern double PDF_q_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);

template <typename Container>
extern double CDF_q_Electron(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);
template <typename Container>
extern double CDF_q_Nucleus(double q, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);

template <typename Container>
extern double Sample_Cos_Theta_Electron(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);
template <typename Container>
extern double Sample_Cos_Theta_Nucleus(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);

template <typename Container>
extern double Sample_q_Electron(std::mt19937& PRNG, double cos_theta, obscura::DM_Particle& DM, double vDM, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);
template <typename Container>
extern double Sample_q_Nucleus(std::mt19937& PRNG, double cos_theta, obscura::DM_Particle& DM, double vDM, obscura::Isotope& nucleus, double nucleus_density, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);

template <typename Container>
extern std::pair<double, double> Sample_Cos_Theta_q_Electron(std::mt19937& PRNG, obscura::DM_Particle& DM, double vDM, double temperature, double number_density_electrons, Container& nuclei, std::vector<double>& number_densities_nuclei, bool use_medium_effects, double qMin);

}	// namespace DaMaSCUS_SUN
#endif