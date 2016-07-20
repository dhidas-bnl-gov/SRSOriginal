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

// UPDATE: Comments

TSurfacePoints_Rectangle::TSurfacePoints_Rectangle ()
{
  // Default Constructor
  // Do nothing

}




TSurfacePoints_Rectangle::TSurfacePoints_Rectangle (int const NX1, int const NX2, TVector3D const& X0, TVector3D const& X1, TVector3D const& X2, int const Normal)
{
  fStartVector = X0;

  fNormal = Normal;

  fNX1 = NX1;
  fNX2 = NX2;
  fNPoints = (size_t) (fNX1 * fNX2);

  double const WidthX1 = (X1 - X0).Mag();
  double const WidthX2 = (X2 - X0).Mag();

  fX1StepSize = WidthX1 / (NX1 - 1);
  fX2StepSize = WidthX2 / (NX2 - 1);

  fX1Vector = (X1 - X0) / (NX1 - 1);
  fX2Vector = (X2 - X0) / (NX2 - 1);

  fNormalVector = fX1Vector.Cross(fX2Vector).UnitVector();

  if (fNormal == -1) {
    fNormalVector *= -1;
  } else if (fNormal != 0 && fNormal != 1) {
    throw;
  }
}






TSurfacePoints_Rectangle::TSurfacePoints_Rectangle (std::string const& p, int const NX1, int const NX2, double const WidthX1, double const WidthX2, TVector3D const& Rotations, TVector3D const& Translation, int const Normal)
{
  // Constructor

  this->Init(p, NX1, NX2, WidthX1, WidthX2, Rotations, Translation, Normal);

}



TSurfacePoints_Rectangle::~TSurfacePoints_Rectangle ()
{
  // Destructor
}




void TSurfacePoints_Rectangle::Init (int const NX1, int const NX2, TVector3D const& X0, TVector3D const& X1, TVector3D const& X2, int const Normal)
{
  fStartVector = X0;

  fNormal = Normal;

  fNX1 = NX1;
  fNX2 = NX2;
  fNPoints = (size_t) (fNX1 * fNX2);

  double const WidthX1 = (X1 - X0).Mag();
  double const WidthX2 = (X2 - X0).Mag();

  fX1StepSize = WidthX1 / (NX1 - 1);
  fX2StepSize = WidthX2 / (NX2 - 1);

  fX1Vector = (X1 - X0) / (NX1 - 1);
  fX2Vector = (X2 - X0) / (NX2 - 1);

  fNormalVector = fX1Vector.Cross(fX2Vector).UnitVector();

  if (fNormal == -1) {
    fNormalVector *= -1;
  } else if (fNormal != 0 && fNormal != 1) {
    throw;
  }


  return;
}






void TSurfacePoints_Rectangle::Init (std::string const& p, int const NX1, int const NX2, double const WidthX1, double const WidthX2, TVector3D const& Rotations, TVector3D const& Translation, int const Normal)
{
  // Constructor

  // I will accept lower-case
  std::string P = p;
  std::transform(P.begin(), P.end(), P.begin(), ::toupper);

  fNX1 = NX1;
  fNX2 = NX2;

  fNormal = Normal;

  fNPoints = (size_t) (fNX1 * fNX2);

  fX1StepSize = WidthX1 / (fNX1 - 1);
  fX2StepSize = WidthX2 / (fNX2 - 1);

  if (P == "XY") {
    fStartVector.SetXYZ(-WidthX1 / 2., -WidthX2 / 2., 0);

    fX1Vector.SetXYZ(fX1StepSize, 0, 0);
    fX2Vector.SetXYZ(0, fX2StepSize, 0);
  } else if (P == "XZ") {
    fStartVector.SetXYZ(-WidthX1 / 2., 0, -WidthX2 / 2.);

    fX1Vector.SetXYZ(fX1StepSize, 0, 0);
    fX2Vector.SetXYZ(0, 0, fX2StepSize);
  } else if (P == "YZ") {
    fStartVector.SetXYZ(0, -WidthX1 / 2., -WidthX2 / 2.);

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

  if (fNormal == -1) {
    fNormalVector *= -1;
  } else if (fNormal != 0 && fNormal != 1) {
    throw;
  }



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
  return fX1StepSize * ix1 - fX1Vector.Mag() * (fNX1 - 1)/ 2.;
}




double TSurfacePoints_Rectangle::GetX2 (size_t const i) const
{
  int const ix2 = i % fNX2;
  return fX2StepSize * ix2 - fX2Vector.Mag() * (fNX2 - 1)/ 2.;
}




double TSurfacePoints_Rectangle::GetElementArea () const
{
  return fX1StepSize * fX2StepSize;
}
