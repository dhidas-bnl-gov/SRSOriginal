#ifndef GUARD_OSCARS_Cuda_h
#define GUARD_OSCARS_Cuda_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Aug 17 08:29:54 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "OSCARS.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <string>

#include "TVector3DC.h"
#include "TField3D_1DRegularized.h"
#include "TSpectrumContainer.h"
#include "TParticleA.h"
#include "TSurfacePoints.h"
#include "T3DScalarContainer.h"

extern "C" int  OSCARS_Cuda_GetDeviceCount ();
extern "C" void OSCARS_Cuda_CalculateFluxGPU (TParticleA& Particle, TSurfacePoints const& Surface, double const Energy_eV, T3DScalarContainer& FluxContainer, int const Dimension = 3, double const Weight = 1, std::string const& OutFileName = "");
extern "C" void OSCARS_Cuda_CalculateSpectrumGPU (TParticleA& Particle, TVector3D const& ObservationPoint, TSpectrumContainer& Spectrum, double const Weight = 1);
extern "C" void OSCARS_Cuda_CalculatePowerDensityGPU (TParticleA& Particle, TSurfacePoints const& Surface, T3DScalarContainer& PowerDensityContainer, int const Dimension, bool const Directional, double const Weight, std::string const& OutFileName = "");







#endif
