#ifndef GUARD_TBField3DZRegularized_h
#define GUARD_TBField3DZRegularized_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Feb 17 18:59:38 EST 2016
//
////////////////////////////////////////////////////////////////////

#include <vector>
#include <array>
#include <string>

#include "TBField3DZ.h"

class TBField3DZRegularized
{
  // This class is designed to be a container for a simple magnetic field.
  // The field array has 4 components, Z and Bx, By, Bz.  You can add elements and sort as needed.
  // If you add elements by hand you also need to call the Sort function yourself

  public:
    TBField3DZRegularized ();
    TBField3DZRegularized (std::string const&);
    TBField3DZRegularized (std::string const&, size_t const);
    TBField3DZRegularized (TBField3DZ&);
    TBField3DZRegularized (TBField3DZ&, size_t const);
    ~TBField3DZRegularized ();

    bool ReadFile (std::string const&);
    bool ReadFileRegularized (std::string const&);
    bool SaveAs (std::string const&, std::string const& Comment = "");

    double GetByAtZ (double const&);
    void   SetZNPointsPerMeter (size_t const);

  private:
    // Parameters to remember in the code.  If fZNPointsPerMeter is zero the grid will default
    // to the grid default defined in TBField3DZ
    double fZFirstPoint;
    double fZLastPoint;
    int    fZNPoints;
    int    fZNPointsPerMeter;
    double fZStepSize;

    // BField.  The vector holds equidistant points.  This is for memory saving and fast lookup.
    std::vector<std::array<double, 3> > fBField;

    bool RegularizeField (TBField3DZ&);
};



















#endif
