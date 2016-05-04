#include "TBFieldSquareWave.h"

#include <cmath>
#include <iostream>


TBFieldSquareWave::TBFieldSquareWave (double const& PeriodLength, double const& NPeriods, double const& CenterZ, double const& MaxBy)
{
  fPeriodLength = PeriodLength;
  fNPeriods = NPeriods;
  fCenterZ = CenterZ;
  fMaxBy = MaxBy;

  fZMin = fCenterZ - (fNPeriods) * (fPeriodLength / 2.);
  fZMax = fCenterZ + (fNPeriods) * (fPeriodLength / 2.);

}





TBFieldSquareWave::~TBFieldSquareWave ()
{
}




double TBFieldSquareWave::GetBx (double const& X, double const& Y, double const& Z) const
{
  return 0;
}




double TBFieldSquareWave::GetBy (double const& X, double const& Y, double const& Z) const
{
  if (Z < fZMin || Z > fZMax) {
    return 0;
  }

  if (Z < fZMin + fPeriodLength || Z > fZMax - fPeriodLength) {
    if (Z < fZMin + fPeriodLength / 2. || Z > fZMax - fPeriodLength / 2.) {
      return 0.25 * fMaxBy * (sin(2. * 3.14 * (Z - fCenterZ) / fPeriodLength) > 0 ? 1.0 : -1.0);
    }
    return 0.75 * fMaxBy * (sin(2. * 3.14 * (Z - fCenterZ) / fPeriodLength) > 0 ? 1.0 : -1.0);
  }

  return fMaxBy * (sin(2. * 3.14 * (Z - fCenterZ) / fPeriodLength) > 0 ? 1.0 : -1.0);
}




double TBFieldSquareWave::GetBz (double const& X, double const& Y, double const& Z) const
{
  return 0;
}
