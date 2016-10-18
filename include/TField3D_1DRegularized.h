#ifndef GUARD_TField3D_1DRegularized_h
#define GUARD_TField3D_1DRegularized_h
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

#include "TField.h"
#include "TField3D_1D.h"


class TField3D_1DRegularized : public TField
{
  // This class is designed to be a container for a simple magnetic field.
  // The field array has 4 components, Z and Bx, By, Bz.  You can add elements and sort as needed.
  // If you add elements by hand you also need to call the Sort function yourself

  public:
    TField3D_1DRegularized ();
    TField3D_1DRegularized (std::string const&);
    TField3D_1DRegularized (std::string const&, TVector3D const&, TVector3D const&, std::vector<double> const&);
    TField3D_1DRegularized (std::string const&, size_t const);
    TField3D_1DRegularized (TField3D_1D&);
    TField3D_1DRegularized (TField3D_1D&, size_t const);
    ~TField3D_1DRegularized ();

    bool ReadFile (std::string const&);
    bool ReadFile (std::string const&, TVector3D const&, TVector3D const&, std::vector<double> const&);
    bool ReadFileRegularized (std::string const&);
    bool SaveAs (std::string const&, std::string const& Comment = "");


    double    GetFx (double const, double const, double const) const;
    double    GetFy (double const, double const, double const) const;
    double    GetFz (double const, double const, double const) const;
    TVector3D GetF  (double const, double const, double const) const;
    TVector3D GetF  (TVector3D const&) const;

    double GetFxAtZ (double const) const;
    double GetFyAtZ (double const) const;
    double GetFzAtZ (double const) const;
    void   SetZNPointsPerMeter (size_t const);

  private:
    // Parameters to remember in the code.  If fZNPointsPerMeter is zero the grid will default
    // to the grid default defined in TField3D_1D
    double fZFirstPoint;
    double fZLastPoint;
    size_t fZNPoints;
    size_t fZNPointsPerMeter;
    double fZStepSize;

    // Field.  The vector holds equidistant points.  This is for memory saving and fast lookup.
    std::vector<std::array<double, 3> > fField;

    bool RegularizeField (TField3D_1D&);
};



















#endif
