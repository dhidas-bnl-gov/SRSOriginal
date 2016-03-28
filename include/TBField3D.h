#ifndef GUARD_TBField3D_h
#define GUARD_TBField3D_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Mar 18 17:19:58 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBField.h"

class TBField3D : public TBField
{

  public:
    TBField3D ();
    ~TBField3D ();

    double GetBx (double const, double const, double const);
    double GetBy (double const, double const, double const);
    double GetBz (double const, double const, double const);

};



















#endif

