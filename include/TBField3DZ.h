#ifndef GUARD_TBField3DZ_h
#define GUARD_TBField3DZ_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Feb 17 18:57:38 EST 2016
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"

#include <vector>
#include <array>
#include <string>

class TBField3DZ
{
  // This class is designed to be a container for a simple magnetic field.
  // The field array has 4 components, Z and Bx, By, Bz.  You can add elements and sort as needed.
  // If you add elements by hand you also need to call the Sort function yourself

  public:
    TBField3DZ ();
    TBField3DZ (std::string const&, TVector3D const& Rotation = TVector3D(0, 0, 0), TVector3D const& Translation = TVector3D(0, 0, 0), std::vector<double> const& Scaling = std::vector<double>());
    ~TBField3DZ ();

    bool Add (double const, double const, double const, double const);
    bool Sort ();
    bool IsSorted ();

    double GetFirstZ ();
    double GetLastZ ();

    bool ReadFile (std::string const&, TVector3D const& Rotation = TVector3D(0, 0, 0), TVector3D const& Translation = TVector3D(0, 0, 0), std::vector<double> const& Scaling = std::vector<double>());
    bool SaveAs (std::string const&, std::string const& Comment = "");

    bool Regularize (std::vector<std::array<double, 3> >&, double&, double&, double&, size_t const NPointsPerMeter = 10000);
    static bool CompareBField3DZ (std::array<double, 4> const&, std::array<double, 4> const&);

  private:
    // BField.  The format is [0] is Z position in Meters, [1, 2, 3] are [Bx, By, Bz] and in Tesla
    std::vector< std::array<double, 4> > fBField;
    bool fIsSorted;
};



















#endif
