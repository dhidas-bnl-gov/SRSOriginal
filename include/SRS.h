#ifndef GUARD_SRS_h
#define GUARD_SRS_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed May 18 20:34:55 EDT 2016
//
// Class to contain all elements and functions for radiation
// simulation.  THIS is the c++ API
//
////////////////////////////////////////////////////////////////////

#include "TSRS.h"

#include <string>

#include "TBFieldContainer.h"
#include "TParticleBeamContainer.h"
#include "TSurfacePoints.h"
#include "TSpectrumContainer.h"


class SRS
{
  // This class is meant to be the main interface to the simulation,
  // also, from all extensions

  public:
    SRS ();
    ~SRS ();


    // Functions related to the magnetic field
    void AddMagneticField (std::string const, std::string const, double const X0 = 0, double const Y0 = 0, double const Z0 = 0);
    double GetBx (double const, double const, double const) const;
    double GetBy (double const, double const, double const) const;
    double GetBz (double const, double const, double const) const;


    // Functions related to the particle beam(s)
    void AddParticleBeam (std::string const&, std::string const&, double const, double const, double const, double const, double const, double const, double const, double const, double const, double const);
    TParticleBeam& GetParticleBeam (std::string const&);
    TParticleA GetNewParticle ();


    // Functions related to Trajectory
    void CalculateTrajectory (TParticleA&);
    TParticleTrajectoryPoints const& GetTrajectory (TParticleA const&) const;

    void SetNPointsTrajectory (size_t const);
    void SetCTStart (double const);
    void SetCTStop  (double const);

    size_t GetNPointsTrajectory () const;
    double GetCTStart () const;
    double GetCTStop  () const;

    void CalculateSpectrum (TParticleA&, TVector3D const&, double const, double const, size_t const, std::string const& OutFileName = "");

    // Power Density calculation
    void CalculatePowerDensity (TParticleA&, TSurfacePoints const&);

    // Flux Calculations
    void CalculateFlux (TParticleA&, TSurfacePoints const&, double const);

  private:
    TBFieldContainer fBFieldContainer;

    TParticleBeamContainer fParticleBeamContainer;

    void Derivatives (double t, double x[], double dxdt[], TParticleA const&);
    void RK4 (double y[], double dydx[], int n, double x, double h, double yout[], TParticleA const&);


    double fCTStart;
    double fCTStop;
    size_t fNPointsTrajectory;


    // Current particle for calculations and rel parameters
    TParticleA fParticle;
    double fCurrent;

    // Spectrum container
    TSpectrumContainer fSpectrum;




};








#endif
