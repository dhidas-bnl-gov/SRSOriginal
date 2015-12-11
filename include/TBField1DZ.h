#ifndef GUARD_TBField1DZ_h
#define GUARD_TBField1DZ_h
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

class TBField1DZ
{
  // This class is designed to be a container for a simple magnetic field.
  // The field array has 2 components, Z and By.  You can add elements and sort as needed.
  // If you add elements by hand you also need to call the Sort function yourself

  public:
    TBField1DZ ();
    TBField1DZ (std::string const&);
    ~TBField1DZ ();

    bool Add (double const, double const);
    bool Sort ();

    bool ReadFile (std::string const&);
    bool SaveAs (std::string const&, std::string const& Comment = "");

    static bool CompareBField1DZ (std::array<double, 2> const&, std::array<double, 2> const&);

  private:
    // BField.  The format is [0] is Z position in Meters, [1] is By in Tesla
    std::vector< std::array<double, 2> > fBField;
};



















#endif
