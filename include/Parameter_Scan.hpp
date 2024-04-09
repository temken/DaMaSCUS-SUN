#ifndef __Parameter_Scan_hpp__
#define __Parameter_Scan_hpp__

#include <vector>

#include "obscura/Configuration.hpp"
#include "obscura/DM_Particle.hpp"
#include "obscura/Direct_Detection.hpp"

#include "Solar_Model.hpp"

namespace DaMaSCUS_SUN
{

// 1. Configuration class for input file, which extends the obscura::Configuration class.

class Configuration : public obscura::Configuration
{
  protected:
	void Import_Parameter_Scan_Parameter();

	virtual void Construct_DM_Particle() override;
	void Construct_DM_Particle_Dark_Photon();

  public:
	std::string run_mode;
	unsigned int isoreflection_rings, interpolation_points;
	unsigned int sample_size, cross_sections;
	double cross_section_min, cross_section_max;
	bool compute_halo_constraints, perform_full_scan, use_medium_effects;
	double zeta;
	explicit Configuration(std::string cfg_filename, int MPI_rank = 0);

	void Print_Summary(int mpi_rank = 0) override;
};

// 2. 	Class to perform parameter scans in the (m_DM, sigma)-plane to search for equal-p-value contours.
//		Either a full scan, or more efficiently and targeted via the square tracing algorithm (STA).

double Compute_p_Value(unsigned int sample_size, obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, unsigned int rate_interpolation_points = 1000, int mpi_rank = 0);

class Parameter_Scan
{
  private:
	std::string results_path;
	std::vector<double> DM_masses;
	std::vector<double> couplings;
	unsigned int sample_size, scattering_rate_interpolation_points;
	double certainty_level;
	std::vector<std::vector<double>> p_value_grid;
	// Check for progress of a previous, incomplete parameter scan to import and continue
	void Import_P_Values();

	// Square tracing algorithm (STA) functions
	bool STA_Point_On_Grid(int row, int column);
	void STA_Go_Forward(int& row, int& column, std::string& STA_direction);
	void STA_Go_Left(int& row, int& column, std::string& STA_direction);
	void STA_Go_Right(int& row, int& column, std::string& STA_direction);
	void STA_Fill_Gaps();

	std::vector<double> Find_Contour_Point(int row, int column, int row_previous, int column_previous, double p_critical);

  public:
	Parameter_Scan(const std::vector<double>& masses, const std::vector<double>& coupl, std::string ID, unsigned int samplesize, unsigned int interpolation_points = 1000, double CL = 0.95);
	Parameter_Scan(Configuration& config);

	void Perform_Full_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank = 0);
	void Perform_STA_Scan(obscura::DM_Particle& DM, obscura::DM_Detector& detector, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, int mpi_rank = 0);

	// Compute the excluded contours based on the p_value_grid using STA to find the contour shape.
	std::vector<std::vector<double>> Limit_Curve();

	void Print_Grid(int mpi_rank = 0, int marker_row = -1, int marker_column = -1);

	void Export_Results(int mpi_rank = 0);
};

}	// namespace DaMaSCUS_SUN
#endif