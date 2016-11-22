////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri May 13 14:52:05 EDT 2016
//
// This is a base class for any surface element.
//
////////////////////////////////////////////////////////////////////

#include "TSurfaceElement_Rectangle.h"

TSurfaceElement_Rectangle::TSurfaceElement_Rectangle ()
{
  // Constructor
}



TSurfaceElement_Rectangle::TSurfaceElement_Rectangle (TSurfacePoint const& P, double const A)
{
  this->Init(P, A);
}




TSurfaceElement_Rectangle::TSurfaceElement_Rectangle (TVector3D const& X, TVector3D const& N, double const A)
{
  this->Init(TSurfacePoint(X, N), A);
}



TSurfaceElement_Rectangle::TSurfaceElement_Rectangle (double const X, double const Y, double const Z, double const XN, double const YN, double const ZN, double const A)
{
  this->Init(TSurfacePoint(X, Y, Z, XN, YN, ZN), A);
}



TSurfaceElement_Rectangle::~TSurfaceElement_Rectangle ()
{
  // Destruction
}




void TSurfaceElement_Rectangle::Init (TSurfacePoint const& P, double const A)
{
  fSurfacePoint = P;
  fArea = A;

  return;
}



TSurfacePoint const& TSurfaceElement_Rectangle::GetSurfacePoint () const
{
  return fSurfacePoint;
}



TVector3D const& TSurfaceElement_Rectangle::GetPoint () const
{
  return fSurfacePoint.GetPoint();
}



TVector3D const& TSurfaceElement_Rectangle::GetNormal () const
{
  return fSurfacePoint.GetNormal();
}



double TSurfaceElement_Rectangle::GetArea () const
{
  return fArea;
}
