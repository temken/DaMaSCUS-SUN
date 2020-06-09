#ifndef __Solar_Model_hpp_
#define __Solar_Model_hpp_

// Headers from libphysica
#include "Linear_Algebra.hpp"
#include "Numerics.hpp"

// Headers from obscura
#include "DM_Particle.hpp"
#include "Target_Nucleus.hpp"

// 1. Nuclear targets in the Sun
class Solar_Nucleus : public obscura::Element
{
  private:
	libphysica::Interpolation number_density;

  public:
	Solar_Nucleus(const obscura::Isotope& isotope, const std::vector<std::vector<double>>& density_table, std::string label = "");
	Solar_Nucleus(const std::vector<obscura::Isotope>& isotopes, const std::vector<std::vector<double>>& density_table, std::string label = "");

	double Number_Density(double r);
	double DM_Scattering_Rate(double r, obscura::DM_Particle& DM, double DM_speed);
};

// 2. Solar model
class Solar_Model
{
  private:
	libphysica::Interpolation mass, temperature, local_escape_speed_squared;
	// Auxiliary functions for the data import
	std::vector<std::vector<double>> raw_data;
	void Import_Raw_Data();
	std::vector<std::vector<double>> Create_Interpolation_Table(unsigned int row) const;
	std::vector<std::vector<double>> Create_Escape_Speed_Table();
	std::vector<std::vector<double>> Create_Number_Density_Table(unsigned int target, double mass) const;
	std::vector<std::vector<double>> Create_Number_Density_Table_Electron();

	// Solar electrons
	libphysica::Interpolation number_density_electron;
	double DM_Scattering_Rate_Electron(double r, obscura::DM_Particle& DM, double DM_speed);

  public:
	std::string name;
	std::vector<Solar_Nucleus> nuclear_targets;

	Solar_Model();

	double Mass(double r);
	double Temperature(double r);
	double Local_Escape_Speed(double r);

	double Number_Density_Nucleus(double r, unsigned int nucleus_index);
	double Number_Density_Electron(double r);

	double Total_DM_Scattering_Rate(double r, obscura::DM_Particle& DM, double DM_speed);

	void Print_Summary(int MPI_rank = 0) const;
};

#endif
