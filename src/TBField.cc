////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Mar 29 08:32:36 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBField.h"


bool TBField::IsWithinRange (double const& X, double const& Y, double const& Z) const
{
  // Is this point in space within the range of the magnetic field definition?

  if (X < fZMin || X > fZMax || X < fYMin || X > fYMax || X < fXMin || X > fXMax) {
    return false;
  }

  return true;
}
