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
#include "T3DScalarContainer.h"


class SRS
{
  // This class is meant to be the main interface to the simulation,
  // also, from all extensions

  public:
    SRS ();
    ~SRS ();


    // Functions related to the magnetic field
    void AddMagneticField (std::string const, std::string const, TVector3D const& R = TVector3D(0, 0, 0), TVector3D const& D = TVector3D(0, 0, 0), std::vector<double> const& S = std::vector<double>());
    void AddMagneticField (TBField*);
    void ClearMagneticFields ();

    double    GetBx (double const, double const, double const) const;
    double    GetBy (double const, double const, double const) const;
    double    GetBz (double const, double const, double const) const;
    TVector3D GetB  (double const, double const, double const) const;
    TVector3D GetB  (TVector3D const&) const;


    // Functions related to the particle beam(s)
    void AddParticleBeam (std::string const&, std::string const&, TVector3D const&, TVector3D const&, double const, double const, double const, double const, double const Charge = 0, double const Mass = 0);
    TParticleBeam& GetParticleBeam (std::string const&);
    size_t GetNParticleBeams () const;
    TParticleA GetNewParticle ();
    TParticleA const&  GetCurrentParticle () const;
    void SetNewParticle ();
    void SetNewParticle (std::string const&, std::string const&);
    void ClearParticleBeams ();


    // Functions related to Trajectory
    void CalculateTrajectory ();
    void CalculateTrajectory (TParticleA&);
    TParticleTrajectoryPoints const& GetTrajectory ();

    void SetNPointsTrajectory (size_t const);
    void SetCTStart (double const);
    void SetCTStop  (double const);
    void SetCTStartStop (double const, double const);

    size_t GetNPointsTrajectory () const;
    double GetCTStart () const;
    double GetCTStop  () const;

    void CalculateSpectrum ();
    void CalculateSpectrum (TVector3D const&, TSpectrumContainer&, double const Weight = 1);
    void CalculateSpectrum (TParticleA&, TVector3D const&, TSpectrumContainer&, double const Weight = 1);
    void CalculateSpectrum (TParticleA&, TVector3D const&, double const, double const, size_t const, std::string const& OutFileName = "");
    void CalculateSpectrum (TVector3D const&, double const, double const, size_t const);
    void CalculateSpectrum (TVector3D const&, std::vector<double> const&);

    TSpectrumContainer const& GetSpectrum () const;

    // Power Density calculation
    void CalculatePowerDensity (TParticleA&, TSurfacePoints const&, int const Dimension = 3, bool const Directional = true, std::string const& OutFileName = "");
    void CalculatePowerDensity (TParticleA&, TSurfacePoints const&, T3DScalarContainer&, int const Dimension = 3, bool const Directional = true, std::string const& OutFileName = "");
    void CalculatePowerDensity (TSurfacePoints const&, T3DScalarContainer&, int const Dimension = 3, bool const Directional = true, std::string const& OutFileName = "");
    double CalculateTotalPower ();
    double CalculateTotalPower (TParticleA&);

    // Flux Calculations
    void CalculateFlux (TParticleA&, TSurfacePoints const&, double const, std::string const& OutFileName = "");
    void CalculateFlux (TParticleA&, TSurfacePoints const&, double const, T3DScalarContainer&, std::string const& OutFileName = "");
    void CalculateFlux (TSurfacePoints const&, double const, T3DScalarContainer&, std::string const& OutFileName = "");

    // Electric Field Calculations
    void CalculateElectricFieldTimeDomain (TVector3D const& Observer, T3DScalarContainer&);
    void CalculateElectricFieldTimeDomain (TVector3D const& Observer, T3DScalarContainer&, TParticleA& Particle);

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
