#ifndef __Solar_Model_hpp_
#define __Solar_Model_hpp_

#include <complex>

#include "libphysica/Linear_Algebra.hpp"
#include "libphysica/Numerics.hpp"

#include "obscura/DM_Particle.hpp"
#include "obscura/Target_Nucleus.hpp"

namespace DaMaSCUS_SUN
{

// 1. Nuclear targets in the Sun
class Solar_Isotope : public obscura::Isotope
{
  private:
	libphysica::Interpolation number_density;

  public:
	Solar_Isotope(const obscura::Isotope& isotope, const std::vector<std::vector<double>>& density_table, double abundance = 1.0);

	double Number_Density(double r);
};

// 3. Solar model
class Solar_Model
{
  private:
	libphysica::Interpolation mass, temperature, local_escape_speed_squared, mass_density;

	// Auxiliary functions for the data import
	std::vector<std::vector<double>> raw_data;
	void Import_Raw_Data();
	std::vector<std::vector<double>> Create_Interpolation_Table(unsigned int row) const;
	std::vector<std::vector<double>> Create_Escape_Speed_Table();
	std::vector<std::vector<double>> Create_Number_Density_Table(unsigned int target, double mass) const;
	std::vector<std::vector<double>> Create_Number_Density_Table_Electron();

	// Solar electrons
	libphysica::Interpolation number_density_electron;

	// Interpolation of total scattering rate
	double interpolation_mass = -1.0, interpolation_coupling_electron = -1.0, interpolation_coupling_nuclei = -1.0;
	bool using_interpolated_rate;
	libphysica::Interpolation_2D rate_interpolation;

	// Medium effects

  public:
	bool use_medium_effects;
	double zeta;
	std::string name;
	std::vector<Solar_Isotope> all_isotopes;
	std::vector<Solar_Isotope> target_isotopes;

	Solar_Model(bool medium_effects = false, double q_cutoff_parameter = 0.0);

	double Mass(double r);
	double Mass_Density(double r);
	double Temperature(double r);
	double Local_Escape_Speed(double r);

	double Number_Density_Nucleus(double r, unsigned int nucleus_index);
	std::vector<double> Number_Densities_Nuclei(double r);
	double Number_Density_Electron(double r);

	double DM_Scattering_Rate_Electron(obscura::DM_Particle& DM, double r, double vDM);
	double DM_Scattering_Rate_Nucleus(obscura::DM_Particle& DM, double r, double vDM, unsigned int nucleus_index);

	double Total_DM_Scattering_Rate(obscura::DM_Particle& DM, double r, double vDM);
	double Total_DM_Scattering_Rate_Computed(obscura::DM_Particle& DM, double r, double vDM);

	double Total_DM_Scattering_Rate_Interpolated(obscura::DM_Particle& DM, double r, double vDM);
	void Interpolate_Total_DM_Scattering_Rate(obscura::DM_Particle& DM, unsigned int N_radius, unsigned int N_speed);

	void Print_Summary(int mpi_rank = 0) const;
};

}	// namespace DaMaSCUS_SUN
#endif
