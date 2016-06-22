#ifndef GUARD_T3DScalarContainer_h
#define GUARD_T3DScalarContainer_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jun 21 17:28:24 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"

#include <vector>

class T3DScalar
{
  public:
    T3DScalar (TVector3D const& X, double const V)
    {
      fX = X;
      fV = V;
    }

    ~T3DScalar ()
    {
    }

    TVector3D const& GetX () const
    {
      return fX;
    }

    double GetV () const
    {
      return fV;
    }

  private:
    TVector3D fX;
    double    fV;
};





class T3DScalarContainer
{
  public:
    T3DScalarContainer ();
    ~T3DScalarContainer ();

    void AddPoint (TVector3D const&, double const);

    size_t GetNPoints () const;

    T3DScalar const& GetPoint (size_t const) const;

  private:
    std::vector<T3DScalar> fValues;
};











#endif
