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

	void Perform_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, int mpi_rank = 0);

	std::vector<std::vector<double>> Limit_Curve(double certainty_level);

	void Import_P_Values(const std::string& file_path);
	void Export_P_Values(const std::string& folder_path, int mpi_rank = 0);
	void Export_Limits(const std::string& folder_path, int mpi_rank = 0, std::vector<double> certainty_levels = {0.9, 0.95, 0.99});
};

#endif