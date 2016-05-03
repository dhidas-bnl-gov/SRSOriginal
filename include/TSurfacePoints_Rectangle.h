#ifndef GUARD_TSurfacePoints_Rectangle_h
#define GUARD_TSurfacePoints_Rectangle_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon May  2 17:55:52 EDT 2016
//
// Surface of a simple rectangle
//
////////////////////////////////////////////////////////////////////

#include "TSurfacePoints.h"

class TSurfacePoints_Rectangle : public TSurfacePoints
{
  public:
    TSurfacePoints_Rectangle (int const&, int const&, double const&, double const&, double const&, double const&, double const&);
    ~TSurfacePoints_Rectangle ();

    TSurfacePoint const GetPoint (size_t const) const;
    size_t GetNPoints () const;

  private:
    int fNX1;
    int fNX2;

    double fX1;
    double fX2;

    double fX0;
    double fY0;
    double fZ0;


};




#endif
