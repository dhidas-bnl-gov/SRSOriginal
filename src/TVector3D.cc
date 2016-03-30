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
  // Probably most used and useful constructor
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
  // Return the X-component
  return fX;
}




double TVector3D::GetY () const
{
  // Return the Y-component
  return fY;
}




double TVector3D::GetZ () const
{
  // Return the Z-component
  return fZ;
}




void TVector3D::SetX (double const& X)
{
  // Set the X component
  fX = X;
  return;
}




void TVector3D::SetY (double const& Y)
{
  // Set the Y component
  fY = Y;
  return;
}




void TVector3D::SetZ (double const& Z)
{
  // Set the Z component
  fZ = Z;
  return;
}




void TVector3D::SetXYZ (double const& X, double const& Y, double const& Z)
{
  // Set the X, Y, and Z components
  fX = X;
  fY = Y;
  fZ = Z;

  return;
}




double TVector3D::Mag() const
{
  // Get the magnitude
  return sqrt(Mag2());
}




double TVector3D::Mag2() const
{
  // Get the magnitude squared
  return fX * fX + fY * fY + fZ * fZ;
}




double TVector3D::Dot(TVector3D const& V) const
{
  // Get the dot product of this dot V
  return fX * V.GetX() + fY * V.GetY() + fZ * V.GetZ();
}




TVector3D TVector3D::Cross (TVector3D const& V) const
{
  // Get the cross product of this cross V using the right hand convention
  return TVector3D(fY * V.GetZ() - V.GetY() * fZ, fZ * V.GetX() - V.GetZ() * fX, fX * V.GetY() - V.GetX() * fY);
}




TVector3D TVector3D::UnitVector () const
{
  // Get a unit vector in the direction of this
  return TVector3D(fX / Mag(), fY / Mag(), fZ / Mag());
}









TVector3D TVector3D::operator + (TVector3D const& V) const
{
  // Vector addition, add components and return a vector
  return TVector3D(fX + V.GetX(), fY + V.GetY(), fZ + V.GetZ());
}




TVector3D TVector3D::operator - (TVector3D const& V) const
{
  // Vector subtraction, subtract components and return a vector
  return TVector3D(fX - V.GetX(), fY - V.GetY(), fZ - V.GetZ());
}




TVector3D TVector3D::operator / (double const& V) const
{
  // Divide vector by some scalar
  return TVector3D(fX / V, fY / V, fZ / V);
}




TVector3D& TVector3D::operator = (TVector3D const& V)
{
  // Assignment operator
  fX = V.GetX();
  fY = V.GetY();
  fZ = V.GetZ();
  return *this;
}




TVector3D TVector3D::operator - ()
{
  // Negative vector
  return TVector3D(-fX, -fY, -fZ);
}




TVector3D& TVector3D::operator += (TVector3D const& V)
{
  // Add a vector to this vector by components
  fX += V.GetX();
  fY += V.GetY();
  fZ += V.GetZ();
  return *this;
}




TVector3D& TVector3D::operator -= (TVector3D const& V)
{
  // Subtract a vector from this vector by components
  fX -= V.GetX();
  fY -= V.GetY();
  fZ -= V.GetZ();
  return *this;
}




TVector3D& TVector3D::operator *= (double const& V)
{
  // Multiply this vector by a scalar
  fX *= V;
  fY *= V;
  fZ *= V;
  return *this;
}




TVector3D& TVector3D::operator /= (double const& V)
{
  // Divide this vector by a scalar
  fX /= V;
  fY /= V;
  fZ /= V;
  return *this;
}




bool TVector3D::operator == (TVector3D const& V) const
{
  // Is this vector equal to V by components
  return fX == V.GetX() && fY == V.GetY() && fZ == V.GetZ();
}




bool TVector3D::operator != (TVector3D const& V) const
{
  // Is any component of this vector not equal to the equivalent component of V
  return fX != V.GetX() || fY != V.GetY() || fZ != V.GetZ();
}









