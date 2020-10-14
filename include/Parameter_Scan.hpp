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

// 1. Configuration class for input file, which extends the obscura::Configuration class.

class Configuration : public obscura::Configuration
{
  protected:
	void Import_Parameter_Scan_Parameter();

  public:
	std::string run_mode;
	unsigned int isoreflection_rings;
	unsigned int sample_size, cross_sections;
	double cross_section_min, cross_section_max;
	bool compute_halo_constraints;
	explicit Configuration(std::string cfg_filename, int MPI_rank = 0);

	void Print_Summary(int mpi_rank = 0) override;
};

// 2. 	Class to perform parameter scans in the (m_DM, sigma)-plane to search for equal-p-value contours.
//		Either a full scan, or more efficiently and targeted via the square tracing algorithm (STA).
class Parameter_Scan
{
  private:
	unsigned int counter = 0;
	std::vector<double> DM_masses;
	std::vector<double> couplings;
	unsigned int sample_size;
	double certainty_level;

	double Compute_p_Value(int row, int col, obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank = 0);

	void STA_Go_Forward(int& row, int& col, std::string& direction);
	void STA_Go_Left(int& row, int& col, std::string& direction);
	void STA_Go_Right(int& row, int& col, std::string& direction);
	void STA_Interpolate_Two_Points(int row, int col, int row_previous, int col_previous, double p_critical);
	void STA_Fill_Gaps();
	bool STA_Is_Point_Within_Bounds(int row, int col);

	bool Compute_Limit_Curve();

  public:
	std::vector<std::vector<double>> p_value_grid, limit_curve;

	Parameter_Scan(Configuration& config);
	Parameter_Scan(const std::vector<double>& masses, const std::vector<double>& coupl, unsigned int samplesize, double CL);

	void Perform_Full_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank = 0);
	void Perform_STA_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, std::string ID, int mpi_rank = 0);

	void Print_Grid(int mpi_rank = 0, int index_coupling = -1, int index_mass = -1);

	void Import_P_Values(const std::string& ID);
	void Export_Results(const std::string& ID, int mpi_rank = 0);
};

}	// namespace DaMaSCUS_SUN
#endif