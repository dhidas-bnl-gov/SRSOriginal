////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jun 30 08:09:53 EDT 2016
//
// UPDATE: Comments
//
////////////////////////////////////////////////////////////////////

#include "TBField3D_Uniform.h"

#include <cmath>

TBField3D_Uniform::TBField3D_Uniform (double const Bx, double const By, double const Bz)
{
  fBField = TVector3D(Bx, By, Bz);
  fWidth  = TVector3D(0, 0, 0);
  fCenter = TVector3D(0, 0, 0);
}




TBField3D_Uniform::TBField3D_Uniform (TVector3D const& BField, TVector3D const& Width, TVector3D const& Center)
{
  fBField = BField;
  fWidth  = Width;
  fCenter = Center;
}


double TBField3D_Uniform::GetBx (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z)).GetX();
}




double TBField3D_Uniform::GetBy (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z)).GetY();
}




double TBField3D_Uniform::GetBz (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z)).GetZ();
}




TVector3D TBField3D_Uniform::GetB (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z));
}




TVector3D TBField3D_Uniform::GetB (TVector3D const& X) const
{
  TVector3D const RX = X - fCenter;

  if (fWidth.GetX() > 0) {
    if (fabs(RX.GetX()) > fWidth.GetX() / 2.) {
      return TVector3D(0, 0, 0);
    }
  }

  if (fWidth.GetY() > 0) {
    if (fabs(RX.GetY()) > fWidth.GetY() / 2.) {
      return TVector3D(0, 0, 0);
    }
  }

  if (fWidth.GetZ() > 0) {
    if (fabs(RX.GetZ()) > fWidth.GetZ() / 2.) {
      return TVector3D(0, 0, 0);
    }
  }

  return fBField;
}

