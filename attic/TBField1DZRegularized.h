#ifndef GUARD_TBField1DZRegularized_h
#define GUARD_TBField1DZRegularized_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Dec 10 17:00:06 EST 2015
//
////////////////////////////////////////////////////////////////////

#include <vector>
#include <array>
#include <string>

#include "TBField1DZ.h"

class TBField1DZRegularized
{
  // This class is designed to be a container for a simple magnetic field.
  // The field array has 2 components, Z and By.  You can add elements and sort as needed.
  // If you add elements by hand you also need to call the Sort function yourself

  public:
    TBField1DZRegularized ();
    TBField1DZRegularized (std::string const&);
    TBField1DZRegularized (std::string const&, size_t const);
    TBField1DZRegularized (TBField1DZ&);
    TBField1DZRegularized (TBField1DZ&, size_t const);
    ~TBField1DZRegularized ();

    bool ReadFile (std::string const&);
    bool ReadFileRegularized (std::string const&);
    bool SaveAs (std::string const&, std::string const& Comment = "");

    double GetByAtZ (double const);
    void   SetZNPointsPerMeter (size_t const);

  private:
    // Parameters to remember in the code.  If fZNPointsPerMeter is zero the grid will default
    // to the grid default defined in TBField1DZ
    double fZFirstPoint;
    double fZLastPoint;
    int    fZNPoints;
    int    fZNPointsPerMeter;
    double fZStepSize;

    // BField.  The vector holds equidistant points.  This is for memory saving and fast lookup.
    std::vector<double> fBField;

    bool RegularizeField (TBField1DZ&);
};



















#endif
