#ifndef GUARD_TSurfacePoints_BBoxSimple_h
#define GUARD_TSurfacePoints_BBoxSimple_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu May 12 17:03:19 EDT 2016
//
//
////////////////////////////////////////////////////////////////////

#include "TSurfacePoints.h"

#include <map>
#include <string>

#include "TSurfacePoints_RectangleSimple.h"



class TSurfacePoints_BBoxSimple : public TSurfacePoints
{
  public:
    TSurfacePoints_BBoxSimple (double const, double const, double const, double const, double const, double const, size_t const, size_t const, size_t const);
    ~TSurfacePoints_BBoxSimple ();

    TSurfacePoint const GetPoint (size_t const) const;
    size_t GetNPoints () const;

    size_t GetSurfaceNumber (size_t const) const;

    double GetX1 (size_t const) const;
    double GetX2 (size_t const) const;

  private:
    // Size of a step in each coord
    double fXStepSize;
    double fYStepSize;
    double fZStepSize;

    // First point on each coord (note: not edge)
    double fXStart;
    double fYStart;
    double fZStart;

    // Length of each edge of box
    double fLX;
    double fLY;
    double fLZ;

    // Center position of box
    double fX0;
    double fY0;
    double fZ0;

    // Number of points in each direction
    size_t fNX;
    size_t fNY;
    size_t fNZ;

    // Number of points on the surfaces perp to direction indicated
    size_t fNPlaneX;
    size_t fNPlaneY;
    size_t fNPlaneZ;

    // Total number of points
    size_t fNPoints;

    // 6 Surfaces on a box (well, at least in 3D..)
    TSurfacePoints_RectangleSimple fSurfaces[6];
};





#endif
