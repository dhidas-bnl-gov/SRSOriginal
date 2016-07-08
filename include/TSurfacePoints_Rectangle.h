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

#include <string>

class TSurfacePoints_Rectangle : public TSurfacePoints
{
  public:
    TSurfacePoints_Rectangle ();
    TSurfacePoints_Rectangle (std::string const&, int const, int const, double const, double const, TVector3D const&, TVector3D const&);
    ~TSurfacePoints_Rectangle ();

    void Init(std::string const&, int const, int const, double const, double const, TVector3D const&, TVector3D const&);
    TSurfacePoint const GetPoint (size_t const) const;
    TVector3D GetXYZ (size_t const) const;
    size_t GetNPoints () const;

    double GetX1 (size_t const) const;
    double GetX2 (size_t const) const;

    double GetElementArea () const;

    enum eInPlane {
      kInPlane_XY,
      kInPlane_XZ,
      kInPlane_YZ
    };

  private:
    int fNX1;
    int fNX2;

    double fWidthX1;
    double fWidthX2;

    double fX1Start;
    double fX2Start;

    double fX1StepSize;
    double fX2StepSize;

    size_t fNPoints;
    eInPlane fInPlane;

    TVector3D fNormalVector;
    TVector3D fX1Vector;
    TVector3D fX2Vector;
    TVector3D fStartVector;


};




#endif
