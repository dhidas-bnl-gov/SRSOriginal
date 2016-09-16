////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jun 30 08:09:53 EDT 2016
//
// UPDATE: Comments
//
////////////////////////////////////////////////////////////////////

#include "TBField3D_UniformBox.h"

#include <cmath>

TBField3D_UniformBox::TBField3D_UniformBox (double const Bx, double const By, double const Bz)
{
  fBField = TVector3D(Bx, By, Bz);
  fWidth  = TVector3D(0, 0, 0);
  fCenter = TVector3D(0, 0, 0);
  fRotated = TVector3D(0, 0, 0);

  fIgnoreAxisX = true;
  fIgnoreAxisY = true;
  fIgnoreAxisZ = true;
}




TBField3D_UniformBox::TBField3D_UniformBox (TVector3D const& BField, TVector3D const& Width, TVector3D const& Center, TVector3D const& Rotations)
{
  fBField = BField;
  fBField.RotateSelfXYZ(Rotations);

  fWidth  = Width;
  fCenter = Center;
  fRotated = Rotations;

  fIgnoreAxisX = false;
  fIgnoreAxisY = false;
  fIgnoreAxisZ = false;

  if (fWidth.GetX() <= 0) {
    fIgnoreAxisX = true;
  }
  if (fWidth.GetY() <= 0) {
    fIgnoreAxisY = true;
  }
  if (fWidth.GetZ() <= 0) {
    fIgnoreAxisZ = true;
  }
}


TBField3D_UniformBox::~TBField3D_UniformBox ()
{
  // Destruction!
}


double TBField3D_UniformBox::GetBx (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z)).GetX();
}




double TBField3D_UniformBox::GetBy (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z)).GetY();
}




double TBField3D_UniformBox::GetBz (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z)).GetZ();
}




TVector3D TBField3D_UniformBox::GetB (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z));
}




TVector3D TBField3D_UniformBox::GetB (TVector3D const& X) const
{

  // Translate back into box frame
  TVector3D XInBoxCoordinates = X;
  XInBoxCoordinates.RotateSelfXYZ(fRotated);

  TVector3D const RX = XInBoxCoordinates - fCenter;

  if (!fIgnoreAxisX && fabs(RX.GetX()) > fabs(fWidth.GetX() / 2.) || !fIgnoreAxisY && fabs(RX.GetY()) > fabs(fWidth.GetY() / 2.) || !fIgnoreAxisZ && fabs(RX.GetZ()) > fabs(fWidth.GetZ() / 2.)) {
    return TVector3D(0, 0, 0);
  }

  return fBField;
}






