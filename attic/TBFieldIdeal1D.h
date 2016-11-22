#ifndef GUARD_TBFieldIdeal1D_h
#define GUARD_TBFieldIdeal1D_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Mar 30 17:22:28 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBField.h"



class TBFieldIdeal1D : public TBField
{
  public:
    TBFieldIdeal1D (double const, double const, double const, double const);
    ~TBFieldIdeal1D ();

    double    GetBx (double const, double const, double const) const;
    double    GetBy (double const, double const, double const) const;
    double    GetBz (double const, double const, double const) const;

    TVector3D GetB (double const, double const, double const) const;
    TVector3D GetB (TVector3D const&) const;


  private:
    double fPeriodLength;
    int    fNPeriods;
    double fCenterZ;
    double fMaxBy;

    double fZMin;
    double fZMax;


};








#endif
