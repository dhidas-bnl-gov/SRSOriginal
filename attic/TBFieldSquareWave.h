#ifndef GUARD_TBFieldSquareWave_h
#define GUARD_TBFieldSquareWave_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Mar 30 17:22:28 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBField.h"



class TBFieldSquareWave : public TBField
{
  public:
    TBFieldSquareWave (double const, double const, double const, double const);
    ~TBFieldSquareWave ();

    double GetBx (double const, double const, double const) const;
    double GetBy (double const, double const, double const) const;
    double GetBz (double const, double const, double const) const;


  private:
    double fPeriodLength;
    int    fNPeriods;
    double fCenterZ;
    double fMaxBy;

    double fZMin;
    double fZMax;


};








#endif
