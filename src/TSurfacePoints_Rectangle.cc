////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jul  7 17:56:06 EDT 2016
//
// Surface of a rectangle
//
////////////////////////////////////////////////////////////////////

#include "TSurfacePoints_Rectangle.h"

#include <algorithm>


TSurfacePoints_Rectangle::TSurfacePoints_Rectangle ()
{
  // Default Constructor
  // Do nothing

}






TSurfacePoints_Rectangle::TSurfacePoints_Rectangle (std::string const& p, int const NX1, int const NX2, double const WidthX1, double const WidthX2, TVector3D const& Rotations, TVector3D const& Translation)
{
  // Constructor

  this->Init(p, NX1, NX2, WidthX1, WidthX2, Rotations, Translation);

}



TSurfacePoints_Rectangle::~TSurfacePoints_Rectangle ()
{
  // Destructor
}




void TSurfacePoints_Rectangle::Init (std::string const& p, int const NX1, int const NX2, double const WidthX1, double const WidthX2, TVector3D const& Rotations, TVector3D const& Translation)
{
  // Constructor

  // I will accept lower-case
  std::string P = p;
  std::transform(P.begin(), P.end(), P.begin(), ::toupper);

  fNX1 = NX1;
  fNX2 = NX2;
  fWidthX1  = WidthX1;
  fWidthX2  = WidthX2;

  fNPoints = (size_t) (fNX1 * fNX2);

  fX1StepSize = fWidthX1 / (fNX1 - 1);
  fX2StepSize = fWidthX2 / (fNX2 - 1);

  if (P == "XY") {
    fStartVector.SetXYZ(-fWidthX1 / 2., -fWidthX2 / 2., 0);

    fX1Vector.SetXYZ(fX1StepSize, 0, 0);
    fX2Vector.SetXYZ(0, fX2StepSize, 0);
  } else if (P == "XZ") {
    fStartVector.SetXYZ(-fWidthX1 / 2., 0, -fWidthX2 / 2.);

    fX1Vector.SetXYZ(fX1StepSize, 0, 0);
    fX2Vector.SetXYZ(0, 0, fX2StepSize);
  } else if (P == "YZ") {
    fStartVector.SetXYZ(0, -fWidthX1 / 2., -fWidthX2 / 2.);

    fX1Vector.SetXYZ(0, fX1StepSize, 0);
    fX2Vector.SetXYZ(0, 0, fX2StepSize);
  } else {
    throw;
  }

  fStartVector.RotateSelfXYZ(Rotations);
  fStartVector += Translation;

  fX1Vector.RotateSelfXYZ(Rotations);
  fX2Vector.RotateSelfXYZ(Rotations);

  fNormalVector = fX1Vector.Cross(fX2Vector).UnitVector();

  return;
}


TSurfacePoint const TSurfacePoints_Rectangle::GetPoint (size_t const i) const
{
  return TSurfacePoint(this->GetXYZ(i), fNormalVector);
}




TVector3D TSurfacePoints_Rectangle::GetXYZ (size_t const i) const
{
  int const ix1 = i / fNX2;
  int const ix2 = i % fNX2;

  return (fStartVector + ix1 * fX1Vector + ix2 * fX2Vector);
}




size_t TSurfacePoints_Rectangle::GetNPoints () const
{
  return fNPoints;
}




double TSurfacePoints_Rectangle::GetX1 (size_t const i) const
{
  int const ix1 = i / fNX2;
  return fX1StepSize * ix1 - fWidthX1 / 2.;
}




double TSurfacePoints_Rectangle::GetX2 (size_t const i) const
{
  int const ix2 = i % fNX2;
  return fX2StepSize * ix2 - fWidthX2 / 2.;
}




double TSurfacePoints_Rectangle::GetElementArea () const
{
  return fX1StepSize * fX2StepSize;
}
