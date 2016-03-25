////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Mar 25 10:24:13 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"

#include <cmath>

TVector3D::TVector3D ()
{
  // Default constructor
}




TVector3D::TVector3D (double const& X, double const& Y, double const& Z)
{
  // Probably most used constructor
  fX = X;
  fY = Y;
  fZ = Z;
}




TVector3D::~TVector3D ()
{
  // Destroy me
}




double TVector3D::GetX () const
{
  return fX;
}




double TVector3D::GetY () const
{
  return fY;
}




double TVector3D::GetZ () const
{
  return fZ;
}




void TVector3D::SetX (double const& X)
{
  fX = X;
  return;
}




void TVector3D::SetY (double const& Y)
{
  fY = Y;
  return;
}




void TVector3D::SetZ (double const& Z)
{
  fZ = Z;
  return;
}




void TVector3D::SetXYZ (double const& X, double const& Y, double const& Z)
{
  fX = X;
  fY = Y;
  fZ = Z;
  return;
}




double TVector3D::Mag() const
{
  return sqrt(Mag2());
}




double TVector3D::Mag2() const
{
  return fX * fX + fY * fY + fZ * fZ;
}




double TVector3D::Dot(TVector3D& V) const
{
  return fX * V.GetX() + fY * V.GetY() + fZ * V.GetZ();
}




TVector3D TVector3D::Cross (TVector3D& V) const
{
  return TVector3D(fY * V.GetZ() - V.GetY() * fZ, fZ * V.GetX() - V.GetZ() * fX, fX * V.GetY() - V.GetX() * fY);
}




TVector3D TVector3D::UnitVector () const
{
  return TVector3D(fX / Mag(), fY / Mag(), fZ / Mag());
}









TVector3D TVector3D::operator + (TVector3D const& V) const
{
  return TVector3D(fX + V.GetX(), fY + V.GetY(), fZ + V.GetZ());
}




TVector3D TVector3D::operator - (TVector3D const& V) const
{
  return TVector3D(fX - V.GetX(), fY - V.GetY(), fZ - V.GetZ());
}




TVector3D TVector3D::operator * (double const& V) const
{
  return TVector3D(fX * V, fY * V, fZ * V);
}




TVector3D TVector3D::operator / (double const& V) const
{
  return TVector3D(fX / V, fY / V, fZ / V);
}




TVector3D& TVector3D::operator = (TVector3D const& V)
{
  fX = V.GetX();
  fY = V.GetY();
  fZ = V.GetZ();
  return *this;
}




TVector3D TVector3D::operator - ()
{
  return TVector3D(-fX, -fY, -fZ);
}




TVector3D& TVector3D::operator += (double const& V)
{
  fX += V;
  fY += V;
  fZ += V;
  return *this;
}




TVector3D& TVector3D::operator -= (double const& V)
{
  fX -= V;
  fY -= V;
  fZ -= V;
  return *this;
}




TVector3D& TVector3D::operator *= (double const& V)
{
  fX *= V;
  fY *= V;
  fZ *= V;
  return *this;
}




TVector3D& TVector3D::operator /= (double const& V)
{
  fX /= V;
  fY /= V;
  fZ /= V;
  return *this;
}




bool TVector3D::operator == (TVector3D const& V) const
{
  return fX == V.GetX() && fY == V.GetY() && fZ == V.GetZ();
}




bool TVector3D::operator != (TVector3D const& V) const
{
  return fX != V.GetX() || fY != V.GetY() || fZ != V.GetZ();
}









