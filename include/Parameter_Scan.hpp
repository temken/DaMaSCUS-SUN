#ifndef __Parameter_Scan_hpp__
#define __Parameter_Scan_hpp__

#include <vector>

// Headers from obscura
#include "Configuration.hpp"
#include "DM_Particle.hpp"
#include "Direct_Detection.hpp"

#include "Solar_Model.hpp"

namespace DaMaSCUS_SUN
{

class Configuration : public obscura::Configuration
{
  protected:
	void Import_Parameter_Scan_Parameter();

  public:
	unsigned int sample_size, cross_sections;
	double cross_section_min, cross_section_max;
	bool compute_halo_constraints;
	explicit Configuration(std::string cfg_filename, int MPI_rank = 0);

	void Print_Summary(int mpi_rank = 0) override;
};

class Parameter_Scan
{
  private:
	std::vector<double> DM_masses;
	std::vector<double> couplings;
	std::vector<std::vector<double>> p_value_grid;
	unsigned int sample_size;

  public:
	Parameter_Scan(Configuration& config);
	Parameter_Scan(const std::vector<double>& masses, const std::vector<double>& coupl, unsigned int samplesize);

	void Perform_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank = 0);

	std::vector<std::vector<double>> Limit_Curve(double certainty_level);

	void Print_Grid(int mpi_rank = 0, int index_coupling = -1, int index_mass = -1);

	void Import_P_Values(const std::string& ID);
	void Export_P_Values(const std::string& ID, int mpi_rank = 0);
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

}	// namespace DaMaSCUS_SUN
#endif