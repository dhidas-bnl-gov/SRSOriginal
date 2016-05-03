#ifndef GUARD_TSurfacePoints_RectangleSimple_h
#define GUARD_TSurfacePoints_RectangleSimple_h
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

#include <string>

class TSurfacePoints_RectangleSimple : public TSurfacePoints
{
  public:
    TSurfacePoints_RectangleSimple (std::string const&, int const&, int const&, double const&, double const&, double const&, double const&, double const&, int const&);
    ~TSurfacePoints_RectangleSimple ();

    TSurfacePoint const GetPoint (size_t const) const;
    size_t GetNPoints () const;

    double GetElementArea () const;

    enum eInPlane {
      kInPlane_XY,
      kInPlane_XZ,
      kInPlane_YZ
    };

  private:
    int fNX1;
    int fNX2;

    double fX1;
    double fX2;

    double fX0;
    double fY0;
    double fZ0;

    double fX1Start;
    double fX2Start;

    double fX1StepSize;
    double fX2StepSize;

    size_t fNPoints;
    eInPlane fInPlane;

    int fNormal;


};




#endif
