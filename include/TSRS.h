#ifndef GUARD_TSRS_h
#define GUARD_TSRS_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed May 18 20:34:55 EDT 2016
//
// Class to contain all elements and functions for radiation
// simulation.
//
// Also: Namespace constants for SRS
//
////////////////////////////////////////////////////////////////////


#include <string>

#include "TBFieldContainer.h"
#include "TParticleBeamContainer.h"
#include "TSurfacePoints.h"


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

    void CalculateSpectrum (TParticleA&, TVector3D const&, double const, double const, size_t const);

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
    // 





};






namespace TSRS {

   /* ************************* */
   /* * Fundamental constants * */
   /* ************************* */

   inline double Pi()       { return 3.14159265358979323846; }
   inline double Pi2()      { return Pi() * Pi(); }
   inline double TwoPi()    { return 2.0 * Pi(); }
   inline double FourPi()   { return 4.0 * Pi(); }
   inline double PiOver2()  { return Pi() / 2.0; }
   inline double PiOver4()  { return Pi() / 4.0; }
   inline double InvPi()    { return 1.0 / Pi(); }
   inline double RadToDeg() { return 180.0 / Pi(); }
   inline double DegToRad() { return Pi() / 180.0; }
   inline double Sqrt2()    { return 1.4142135623730950488016887242097; }
   inline double Sqrt2Pi()  { return 2.5066282746310002416123552393401; }

   // e (base of natural log)
   inline double E()        { return 2.71828182845904523536; }

   // base-10 log of e  (to convert ln to log)
   inline double LogE()     { return 0.43429448190325182765; }

   // velocity of light
   inline double C()        { return 2.99792458e8; }        // m s^-1

   // Planck's constant
   inline double H()        { return 6.62606876e-34; }      // J s

   // h-bar (h over 2 pi)
   inline double Hbar()     { return 1.054571596e-34; }     // J s

   // hc (h * c)
   inline double HC()       { return H() * C(); }           // J m



   // Elementary charge
   inline double Qe()       { return 1.602176462e-19; }     // C
   inline double Me()       { return 9.10938356e-31; }      // kg


   // Elementary charge over mass of electron
   inline double QeOverMe() { return Qe() / Me(); }         // C kg^-1

   // Permitivity of vacuum
   inline double Epsilon0() { return 8.854187817E-12; }     // Add units
   inline double Mu0()      { return 1.2566370614E-6; }     // Add units

   inline double FrequencyToEv         (double const f) { return f * H() / Qe();    } // eV
   inline double AngularFrequencyToEv  (double const w) { return w * Hbar() / Qe(); } // eV
   inline double EvToAngularFrequency  (double const e) { return e * Qe() / Hbar(); } // rad s^-1
   inline double EvToFrequency         (double const e) { return e * Qe() / H();    } // s^-1
   inline double kgToGeV               (double const m) { return 1e-9 * m * C() * C() / Qe(); } // GeV 

}





#endif
