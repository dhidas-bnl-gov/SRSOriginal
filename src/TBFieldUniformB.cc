////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon Mar 28 10:56:20 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBFieldUniformB.h"


TBFieldUniformB::TBFieldUniformB (double const& X, double const& Y, double const& Z)
{
  // Constructor
  fBx = X;
  fBy = Y;
  fBz = Z;
}




TBFieldUniformB::~TBFieldUniformB ()
{
  // Destructor
}




double TBFieldUniformB::GetBx (double const& X, double const& Y, double const& Z) const
{
  // Get the X component of the field
  return fBx;
}




double TBFieldUniformB::GetBy (double const& X, double const& Y, double const& Z) const
{
  // Get the Y component of the field
  return fBy;
}




double TBFieldUniformB::GetBz (double const& X, double const& Y, double const& Z) const
{
  // Get the Z component of the field
  return fBz;
}
