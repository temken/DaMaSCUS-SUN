#ifndef __Reflection_Spectrum_hpp__
#define __Reflection_Spectrum_hpp__

#include "libphysica/Numerics.hpp"

#include "obscura/DM_Distribution.hpp"

#include "Data_Generation.hpp"
#include "Solar_Model.hpp"

namespace DaMaSCUS_SUN
{

class Reflection_Spectrum : public obscura::DM_Distribution
{
  private:
	libphysica::Interpolation kde_speed;
	double total_entering_rate, total_reflection_rate, total_reflection_flux;
	double distance;

  public:
	//Constructors
	Reflection_Spectrum(const Simulation_Data& simulation_data, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, double mDM, int iso_ring = 0);

	virtual double PDF_Speed(double v) override;

	double Differential_Spectrum(double v);
	virtual double Differential_DM_Flux(double v, double mDM = 0.0) override;

	void Set_Distance(double d);

	virtual void Print_Summary(int mpi_rank = 0) override;
};

double DM_Entering_Rate(Solar_Model& solar_model, obscura::DM_Distribution& halo_model, double mDM);

}	// namespace DaMaSCUS_SUN

#endif