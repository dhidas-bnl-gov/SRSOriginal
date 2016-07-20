////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jul 19 08:13:02 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBField3D_IdealUndulator.h"

#include "TSRS.h"
#include <cmath>


TBField3D_IdealUndulator::TBField3D_IdealUndulator ()
{
  // Constructor
}




TBField3D_IdealUndulator::TBField3D_IdealUndulator (TVector3D const& BField, TVector3D const& Period, int const NPeriods, TVector3D const& Center, double const Phase)
{
  // Typical constructor
  this->Init(BField, Period, NPeriods, Center, Phase);
}




TBField3D_IdealUndulator::~TBField3D_IdealUndulator ()
{
  // Destruction!
}




void TBField3D_IdealUndulator::Init (TVector3D const& BField, TVector3D const& Period, int const NPeriods, TVector3D const& Center, double const Phase)
{
  fBField   = BField;
  fPeriod   = Period;
  fNPeriods = NPeriods;
  fCenter   = Center;
  fPhase    = Phase;

  fPeriodLength = fPeriod.Mag();
  fPeriodUnitVector = fPeriod.UnitVector();

  fUndulatorLength = fPeriod.Mag() * (fNPeriods + 2);

  std::cout << "BField " << fBField << std::endl;
  std::cout << "Period " << fPeriod << std::endl;
  std::cout << "NPeriods " << fNPeriods << std::endl;
  std::cout << "Center " << fCenter << std::endl;
  std::cout << "Phase " << fPhase << std::endl;
  std::cout << "UndulatorLength " << fUndulatorLength << std::endl;
  std::cout << "PeriodLength " << fPeriodLength << std::endl;
  std::cout << "PeriodUnitVector " << fPeriodUnitVector << std::endl;


  return;
}




double TBField3D_IdealUndulator::GetBx (double const X, double const Y, double const Z) const
{
  return this->GetB(X, Y, Z).GetX();
}




double TBField3D_IdealUndulator::GetBy (double const X, double const Y, double const Z) const
{
  return this->GetB(X, Y, Z).GetY();
}




double TBField3D_IdealUndulator::GetBz (double const X, double const Y, double const Z) const
{
  return this->GetB(X, Y, Z).GetZ();
}




TVector3D TBField3D_IdealUndulator::GetB (TVector3D const& X) const
{

  // UPDATE: Check min.max limits here?  maybe not

  // How far are you from the "center" in the correct direction
  double const D = (X - fCenter).Dot( fPeriodUnitVector );

  //std::cout << "D " << D << std::endl;

  double const PhaseShift = fPhase * fPeriod.Mag() / TSRS::TwoPi ();
  //std::cout << "PhaseShift " << PhaseShift << std::endl;

  TVector3D B(0, 0, 0);


  if (D > fUndulatorLength / 2. + PhaseShift || D < -fUndulatorLength / 2. + PhaseShift) {
    return B;
  }

  if (D < -fUndulatorLength / 2. + PhaseShift + fPeriodLength || D > fUndulatorLength / 2. + PhaseShift - fPeriodLength) {
    if (D < -fUndulatorLength / 2. + PhaseShift + fPeriodLength / 2. || D > fUndulatorLength / 2. + PhaseShift - fPeriodLength / 2.) {

      return 0.25 * fBField * sin(TSRS::TwoPi() * (D - PhaseShift) / fPeriodLength);
    } 

    return 0.75 * fBField * sin(TSRS::TwoPi() * (D - PhaseShift) / fPeriodLength);
    //return 0.75 * fBField * sin(TSRS::TwoPi() * (D + PhaseShift) / fPeriodLength);
  } 

  return fBField * sin(TSRS::TwoPi() * (D - PhaseShift) / fPeriodLength);
  //return fBField * sin(TSRS::TwoPi() * (D + PhaseShift) / fPeriodLength);
}




TVector3D TBField3D_IdealUndulator::GetB (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z));
}


