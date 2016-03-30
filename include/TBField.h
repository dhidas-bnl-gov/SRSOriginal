#ifndef GUARD_TBField_h
#define GUARD_TBField_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Mar 18 17:13:44 EDT 2016
//
////////////////////////////////////////////////////////////////////


class TBField
{
  // This class is designed to be a base class for a magnetic field object.

  public:
    //virtual ~TBField () = 0;

    virtual double GetBx (double const&, double const&, double const&) const = 0;
    virtual double GetBy (double const&, double const&, double const&) const = 0;
    virtual double GetBz (double const&, double const&, double const&) const = 0;

    bool IsWithinRange (double const&, double const&, double const&) const;

  private:
    double fXMin;
    double fXMax;
    double fYMin;
    double fYMax;
    double fZMin;
    double fZMax;

};



















#endif
