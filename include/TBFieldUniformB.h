#ifndef GUARD_TBFieldUniformB_h
#define GUARD_TBFieldUniformB_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon Mar 28 10:56:20 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBField.h"



class TBFieldUniformB : public TBField
{
  public:
    TBFieldUniformB (double const, double const, double const);
    ~TBFieldUniformB ();

    double    GetBx (double const, double const, double const) const;
    double    GetBy (double const, double const, double const) const;
    double    GetBz (double const, double const, double const) const;
    TVector3D GetB  (double const, double const, double const) const;

  private:
    double fBx;
    double fBy;
    double fBz;
};







#endif
