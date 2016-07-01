#ifndef GUARD_TBField3D_Uniform_h
#define GUARD_TBField3D_Uniform_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jun 30 08:09:53 EDT 2016
//
// UPDATE: Comments
//
////////////////////////////////////////////////////////////////////

#include "TBField.h"
#include "TVector3D.h"

class TBField3D_Uniform : public TBField
{
  public:
    TBField3D_Uniform (double const, double const, double const);
    TBField3D_Uniform (TVector3D const&, TVector3D const& Width = TVector3D(0, 0, 0), TVector3D const& Center = TVector3D(0, 0, 0));
    ~TBField3D_Uniform ();

    double    GetBx (double const, double const, double const) const;
    double    GetBy (double const, double const, double const) const;
    double    GetBz (double const, double const, double const) const;
    TVector3D GetB  (double const, double const, double const) const;
    TVector3D GetB  (TVector3D const&) const;

  private:
    TVector3D fBField;
    TVector3D fWidth;
    TVector3D fCenter;
};






#endif

