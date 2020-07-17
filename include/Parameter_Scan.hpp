#ifndef __Parameter_Scan_hpp__
#define __Parameter_Scan_hpp__

#include <vector>

// Headers from obscura
#include "DM_Particle.hpp"
#include "Direct_Detection.hpp"

#include "Solar_Model.hpp"

class Parameter_Scan
{
  private:
	std::vector<double> DM_masses;
	std::vector<double> couplings;
	std::vector<std::vector<double>> p_value_grid;
	unsigned int sample_size;

  public:
	Parameter_Scan(const std::vector<double> masses, const std::vector<double>& coupl, unsigned int samplesize);

	void Perform_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank = 0);

	std::vector<std::vector<double>> Limit_Curve(double certainty_level);

	void Import_P_Values(const std::string& file_path);
	void Export_P_Values(const std::string& folder_path, int mpi_rank = 0);
	void Export_Limits(const std::string& folder_path, int mpi_rank = 0, std::vector<double> certainty_levels = {0.9, 0.95, 0.99});
};

class Solar_Reflection_Limit
{
  private:
	unsigned int sample_size;
	double coupling_min, coupling_max;
	std::vector<double> masses;
	std::vector<double> limits;
	double certainty_level;

  public:
	Solar_Reflection_Limit(unsigned int Nsample, double mMin, double mMax, unsigned int Nmass, double c_min, double c_max, double CL = 0.95);

	double Upper_Limit(double mass, obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank = 0);
	void Compute_Limit_Curve(std::string& ID, obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank = 0);
	void Export_Curve(std::string& ID, int mpi_rank = 0);
};

#endif