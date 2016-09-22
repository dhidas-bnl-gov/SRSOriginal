#ifndef GUARD_TField3DGrid_h
#define GUARD_TField3DGrid_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Sep 20 07:55:00 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TField.h"

#include <string>
#include <vector>

class TField3DGrid : public TField
{

  public:
    TField3DGrid ();
    TField3DGrid (std::string const&, std::string const& FileFormat = "");
    ~TField3DGrid ();

    double GetFx (double const, double const, double const) const;
    double GetFy (double const, double const, double const) const;
    double GetFz (double const, double const, double const) const;
    TVector3D GetF (double const, double const, double const) const;
    TVector3D GetF (TVector3D const&) const;

    size_t GetIndex (size_t const, size_t const, size_t const) const;

    double GetHeaderValueSRW (std::string const&, const char CommentChar = '#') const;

    void ReadFile_SRW (std::string const&);
    void ReadFile_SPECTRA (std::string const&);

    enum TField3DGrid_DIMX {
      kDIMX_X,
      kDIMX_Y,
      kDIMX_Z,
      kDIMX_XY,
      kDIMX_XZ,
      kDIMX_YZ,
      kDIMX_XYZ
    };



  private:
    // Dimension and position data
    int    fNX;
    int    fNY;
    int    fNZ;
    double fXStart;
    double fYStart;
    double fZStart;
    double fXStep;
    double fYStep;
    double fZStep;
    double fXStop;
    double fYStop;
    double fZStop;

    bool fHasX;
    bool fHasY;
    bool fHasZ;
    int  fXDIM;

    TField3DGrid_DIMX fDIMX;

    // Field data
    std::vector<TVector3D> fData;

};



















#endif

