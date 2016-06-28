#ifndef GUARD_TBFieldPythonFunction_h
#define GUARD_TBFieldPythonFunction_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon Jun 27 07:54:27 EDT 2016
//
////////////////////////////////////////////////////////////////////
#include "Python.h"

#include "TBField.h"

class TBFieldPythonFunction : public TBField
{
  public:
    TBFieldPythonFunction (PyObject*);
    ~TBFieldPythonFunction ();

    double    GetBx (double const, double const, double const) const;
    double    GetBy (double const, double const, double const) const;
    double    GetBz (double const, double const, double const) const;
    TVector3D GetB  (double const, double const, double const) const;
    TVector3D GetB  (TVector3D const&) const;

  private:
    PyObject* fPythonFunction;

};







#endif
