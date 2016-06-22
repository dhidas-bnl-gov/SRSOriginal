////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon May  2 17:55:52 EDT 2016
//
// Surface of a simple rectangle
//
////////////////////////////////////////////////////////////////////

#include "TSurfacePoints_RectangleSimple.h"

#include <algorithm>


TSurfacePoints_RectangleSimple::TSurfacePoints_RectangleSimple ()
{
  // Default Constructor
  // Do nothing

}






TSurfacePoints_RectangleSimple::TSurfacePoints_RectangleSimple (std::string const& p, int const nx1, int const nx2, double const x1, double const x2, TVector3D const& Center, int const n)
{
  // Constructor

  this->Init(p, nx1, nx2, x1, x2, Center, n);

}



TSurfacePoints_RectangleSimple::TSurfacePoints_RectangleSimple (std::string const& p, int const nx1, int const nx2, double const x1, double const x2, double const x, double const y, double const z, int const n)
{
  // Constructor

  this->Init(p, nx1, nx2, x1, x2, TVector3D(x, y, z), n);

}



TSurfacePoints_RectangleSimple::~TSurfacePoints_RectangleSimple ()
{
  // Destructor
}




void TSurfacePoints_RectangleSimple::Init (std::string const& p, int const nx1, int const nx2, double const x1, double const x2, double const x, double const y, double const z, int const n)
{
  this->Init(p, nx1, nx2, x1, x2, TVector3D(x, y, z), n);
  return;
}




void TSurfacePoints_RectangleSimple::Init (std::string const& p, int const nx1, int const nx2, double const x1, double const x2, TVector3D const& Center, int const n)
{
  // Constructor

  // I will accept lower-case
  std::string P = p;
  std::transform(P.begin(), P.end(), P.begin(), ::toupper);

  fNX1 = nx1;
  fNX2 = nx2;
  fX1  = x1;
  fX2  = x2;

  fX0 = Center.GetX();
  fY0 = Center.GetY();
  fZ0 = Center.GetZ();

  fNPoints = (size_t) (nx1 * nx2);

  fX1StepSize = fX1 / (fNX1 - 1);
  fX2StepSize = fX2 / (fNX2 - 1);

  if (P == "XY") {
    fInPlane = kInPlane_XY;
    fX1Start = fX0 - fX1 / 2.;
    fX2Start = fY0 - fX2 / 2.;
  } else if (P == "XZ") {
    fInPlane = kInPlane_XZ;
    fX1Start = fX0 - fX1 / 2.;
    fX2Start = fZ0 - fX2 / 2.;
  } else if (P == "YZ") {
    fInPlane = kInPlane_YZ;
    fX1Start = fY0 - fX1 / 2.;
    fX2Start = fZ0 - fX2 / 2.;
  } else {
    throw;
  }

  if (n == 1) {
    fNormal = 1;
  } else if (n == -1) {
    fNormal = -1;
  } else {
    throw;
  }


}


TSurfacePoint const TSurfacePoints_RectangleSimple::GetPoint (size_t const i) const
{
  int const ix1 = i / fNX2;
  int const ix2 = i % fNX2;

  if (fInPlane == kInPlane_XY) {
    return TSurfacePoint(fX1StepSize * ix1 + fX1Start, fX2StepSize * ix2 + fX2Start, fZ0, 0, 0, fNormal);
  } else if (fInPlane == kInPlane_XZ) {
    return TSurfacePoint(fX1StepSize * ix1 + fX1Start, fY0, fX2StepSize * ix2 + fX2Start, 0, fNormal, 0);
  } else if (fInPlane == kInPlane_YZ) {
    return TSurfacePoint(fX0, fX1StepSize * ix1 + fX1Start, fX2StepSize * ix2 + fX2Start, fNormal, 0, 0);
  } else {
    throw;
  }


  TSurfacePoint p;
  return p;
}




size_t TSurfacePoints_RectangleSimple::GetNPoints () const
{
  return fNPoints;
}




double TSurfacePoints_RectangleSimple::GetX1 (size_t const i) const
{
  int const ix1 = i / fNX2;
  return fX1StepSize * ix1 + fX1Start;
}




double TSurfacePoints_RectangleSimple::GetX2 (size_t const i) const
{
  int const ix2 = i % fNX2;
  return fX2StepSize * ix2 + fX2Start;
}




double TSurfacePoints_RectangleSimple::GetElementArea () const
{
  return fX1StepSize * fX2StepSize;
}
