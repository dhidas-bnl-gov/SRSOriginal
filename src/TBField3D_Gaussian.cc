////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Jun 29 08:52:22 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBField3D_Gaussian.h"
#include "TSRS.h"

#include <cmath>

TBField3D_Gaussian::TBField3D_Gaussian ()
{
  // Constructor
}



TBField3D_Gaussian::TBField3D_Gaussian (TVector3D const& BField, TVector3D const& Center, TVector3D const& Sigma)
{
  // Constructor you should use.. just a suggestion...

  fBField = BField;
  fCenter = Center;
  fSigma = Sigma;
}



double TBField3D_Gaussian::GetBx (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z)).GetX();
}




double TBField3D_Gaussian::GetBy (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z)).GetY();
}




double TBField3D_Gaussian::GetBz (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z)).GetZ();
}




TVector3D TBField3D_Gaussian::GetB (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z));
}




TVector3D TBField3D_Gaussian::GetB (TVector3D const& X) const
{
  //Fraction *=  1. / (TSRS::TwoPi() * fSigma.GetX() *fSigma.GetX()) * exp( -pow(X.GetX() - fCenter.GetX(), 2) / (2 * fSigma.GetX() * fSigma.GetX()) );
  double Fraction = 1;
  if (fSigma.GetX() > 0) {
    Fraction *= exp(-pow((X.GetX() - fCenter.GetX()) / fSigma.GetX(), 2) / 2.);
  }
  if (fSigma.GetY() > 0) {
    Fraction *= exp(-pow((X.GetY() - fCenter.GetY()) / fSigma.GetY(), 2) / 2.);
  }
  if (fSigma.GetZ() > 0) {
    Fraction *= exp(-pow((X.GetZ() - fCenter.GetZ()) / fSigma.GetZ(), 2) / 2.);
  }

  return Fraction * fBField;
}
