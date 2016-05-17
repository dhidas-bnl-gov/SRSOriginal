#include "TBFieldIdeal1D.h"

#include <cmath>
#include <iostream>


TBFieldIdeal1D::TBFieldIdeal1D (double const& PeriodLength, double const& NPeriods, double const& CenterZ, double const& MaxBy)
{
  fPeriodLength = PeriodLength;
  fNPeriods = NPeriods;
  fCenterZ = CenterZ;
  fMaxBy = MaxBy;

  fZMin = fCenterZ - (fNPeriods + 1) * (fPeriodLength / 2.);
  fZMax = fCenterZ + (fNPeriods + 1) * (fPeriodLength / 2.);

}





TBFieldIdeal1D::~TBFieldIdeal1D ()
{
}




double TBFieldIdeal1D::GetBx (double const& X, double const& Y, double const& Z) const
{
  return 0;
}




double TBFieldIdeal1D::GetBy (double const& X, double const& Y, double const& Z) const
{
  if (Z < fZMin || Z > fZMax) {
    return 0;
  }


  if (Z < fZMin + fPeriodLength || Z > fZMax - fPeriodLength) {

    if (Z < fZMin + fPeriodLength / 2. || Z > fZMax - fPeriodLength / 2.) {
      return 0.25 * fMaxBy * sin(2. * 3.14 * (Z - fCenterZ) / fPeriodLength);
    }

    return 0.75 * fMaxBy * sin(2. * 3.14 * (Z - fCenterZ) / fPeriodLength);
  }

  return fMaxBy * sin(2. * 3.14 * (Z - fCenterZ) / fPeriodLength);
}




double TBFieldIdeal1D::GetBz (double const& X, double const& Y, double const& Z) const
{
  return 0;
}
