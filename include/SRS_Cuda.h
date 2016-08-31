#ifndef GUARD_SRS_Cuda_h
#define GUARD_SRS_Cuda_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Aug 17 08:29:54 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "SRS.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <string>

#include "TVector3DC.h"
#include "TBField1DZRegularized.h"
#include "TBField3DZRegularized.h"
#include "TSpectrumContainer.h"
#include "TParticleA.h"
#include "TSurfacePoints.h"
#include "T3DScalarContainer.h"

void SRS_Cuda_CalculateFluxGPU (TParticleA& Particle, TSurfacePoints const& Surface, double const Energy_eV, T3DScalarContainer& FluxContainer, int const Dimension = 3, double const Weight = 1, std::string const& OutFileName = "");
void SRS_Cuda_CalculateSpectrumGPU (TParticleA& Particle, TVector3D const& ObservationPoint, TSpectrumContainer& Spectrum, double const Weight = 1);
void SRS_Cuda_CalculatePowerDensityGPU (TParticleA& Particle, TSurfacePoints const& Surface, T3DScalarContainer& PowerDensityContainer, int const Dimension, bool const Directional, std::string const& OutFileName = "");







#endif
