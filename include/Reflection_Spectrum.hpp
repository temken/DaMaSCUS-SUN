#ifndef __Reflection_Spectrum_hpp__
#define __Reflection_Spectrum_hpp__

// Headers from libphysica
#include "Numerics.hpp"

// Headers from obscura
#include "DM_Distribution.hpp"

#include "Data_Generation.hpp"
#include "Solar_Model.hpp"

class Reflection_Spectrum : public obscura::DM_Distribution
{
  private:
	libphysica::Interpolation spectrum_interpolation;

  public:
	//Constructors
	Reflection_Spectrum(const Simulation_Data& simulation_data, unsigned int iso_ring = 0);

	virtual double PDF_Speed(double v) override;

	double Differential_Spectrum(double v);
	double Differential_Flux(double v);

	virtual void Print_Summary(int mpi_rank = 0) override;
};

double DM_Entering_Rate(Solar_Model& solar_model, obscura::DM_Distribution& halo_model, double mDM);
#endif