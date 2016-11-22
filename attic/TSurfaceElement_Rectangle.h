#ifndef GUARD_TSurfaceElement_Rectangle_h
#define GUARD_TSurfaceElement_Rectangle_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri May 13 14:52:05 EDT 2016
//
// This is a base class for any surface element.
//
////////////////////////////////////////////////////////////////////

#include "TSurfaceElement.h"

class TSurfaceElement_Rectangle : public TSurfaceElement
{
  public:

    TSurfaceElement_Rectangle ();
    TSurfaceElement_Rectangle (TSurfacePoint const&, double const);
    TSurfaceElement_Rectangle (TVector3D const&, TVector3D const&, double const);
    TSurfaceElement_Rectangle (double const, double const, double const, double const, double const, double const, double const);
    ~TSurfaceElement_Rectangle ();

    void Init (TSurfacePoint const&, double const);

    TSurfacePoint const& GetSurfacePoint () const;
    TVector3D     const& GetPoint () const ;
    TVector3D     const& GetNormal () const;
    double        GetArea () const;

};







#endif

