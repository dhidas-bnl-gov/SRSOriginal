////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Sep 20 07:55:30 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TField3DGrid.h"

#include <fstream>
#include <sstream>

TField3DGrid::TField3DGrid ()
{
  fRotated.SetXYZ(0, 0, 0);
  fTranslation.SetXYZ(0, 0, 0);
}




TField3DGrid::TField3DGrid (std::string const& InFileName, std::string const& FileFormat, TVector3D const& Rotations, TVector3D const& Translation, char const CommentChar)
{
  // I will accept lower-case
  std::string format = FileFormat;
  std::transform(format.begin(), format.end(), format.begin(), ::toupper);

  // Which file format are you looking at?
  if (format == "OSCARS") {
    this->ReadFile(InFileName, Rotations, Translation);
  } else if (format == "SPECTRA") {
    this->ReadFile_SPECTRA(InFileName, Rotations, Translation, CommentChar);
  } else if (format == "SRW") {
    this->ReadFile_SRW(InFileName, Rotations, Translation, CommentChar);
  } else {
  }
}




TField3DGrid::~TField3DGrid ()
{
}




double TField3DGrid::GetFx (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z)).GetX();
}




double TField3DGrid::GetFy (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z)).GetY();
}




double TField3DGrid::GetFz (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z)).GetZ();
}





TVector3D TField3DGrid::GetF (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z));
}



size_t TField3DGrid::GetIndex (size_t const ix, size_t const iy, size_t const iz) const
{
  return ix * fNY * fNZ + iy * fNZ + iz;
}



TVector3D TField3DGrid::GetF (TVector3D const& XIN) const
{
  // Get the field at a point in space.  Must rotate point into coordinate system, then translate it.

  // Rotate and Translate
  TVector3D X(XIN);
  X.RotateSelfXYZ(fRotated);
  X -= fTranslation;

  // If outside the range, return a zero
  if (fNX > 1 && (X.GetX() < fXStart || X.GetX() > fXStop)) {
    return TVector3D(0, 0, 0);
  }
  if (fNY > 1 && (X.GetY() < fYStart || X.GetY() > fYStop)) {
    return TVector3D(0, 0, 0);
  }
  if (fNZ > 1 && (X.GetZ() < fZStart || X.GetZ() > fZStop)) {
    return TVector3D(0, 0, 0);
  }
      
  // Get index in each dimension relative to start and stop
  size_t const nx = fNX > 1 ? (X.GetX() - fXStart) / fXStep      : 0;
  double const dx = fNX > 1 ? (X.GetX() - fXStart) - nx * fXStep : 0;

  size_t const ny = fNY > 1 ? (X.GetY() - fYStart) / fYStep      : 0;
  double const dy = fNY > 1 ? (X.GetY() - fYStart) - ny * fYStep : 0;

  size_t const nz = fNZ > 1 ? (X.GetZ() - fZStart) / fZStep      : 0;
  double const dz = fNZ > 1 ? (X.GetZ() - fZStart) - nz * fZStep : 0;

  switch (fDIMX) {
    case kDIMX_XYZ:
      {
        // First move in X to find the 4 points of square
        size_t const i000 = GetIndex(nx + 0, ny + 0, nz + 0);
        size_t const i100 = GetIndex(nx + 1, ny + 0, nz + 0);
        TVector3D const v00 = fData[i000] + dx * (fData[i100] - fData[i000]) / fXStep;

        size_t const i010 = GetIndex(nx + 0, ny + 1, nz + 0);
        size_t const i110 = GetIndex(nx + 1, ny + 1, nz + 0);
        TVector3D const v10 = fData[i010] + dx * (fData[i110] - fData[i010]) / fXStep;

        size_t const i001 = GetIndex(nx + 0, ny + 0, nz + 1);
        size_t const i101 = GetIndex(nx + 1, ny + 0, nz + 1);
        TVector3D const v01 = fData[i001] + dx * (fData[i101] - fData[i001]) / fXStep;

        size_t const i011 = GetIndex(nx + 0, ny + 1, nz + 1);
        size_t const i111 = GetIndex(nx + 1, ny + 1, nz + 1);
        TVector3D const v11 = fData[i011] + dx * (fData[i111] - fData[i011]) / fXStep;


        // Step in Y to find 2 points
        TVector3D const v0 = v00 + dy * (v10 - v00) / fYStep;
        TVector3D const v1 = v01 + dy * (v11 - v01) / fYStep;

        // Step in Z to find point
        return v0 + dz * (v1 - v0) / fZStep;
      }
      break;
    case kDIMX_X:
      {
        size_t const i0 = nx + 0;
        size_t const i1 = nx + 1;
        return fData[i0] + dx * (fData[i1] - fData[i0]) / fXStep;
      }
    case kDIMX_Y:
      {
        size_t const i0 = ny + 0;
        size_t const i1 = ny + 1;
        return fData[i0] + dy * (fData[i1] - fData[i0]) / fYStep;
      }
    case kDIMX_Z:
      {
        size_t const i0 = nz + 0;
        size_t const i1 = nz + 1;
        return fData[i0] + dz * (fData[i1] - fData[i0]) / fZStep;
      }
    case kDIMX_XY:
      {
        size_t const i00 = GetIndex(nx + 0, ny + 0, 0);
        size_t const i10 = GetIndex(nx + 1, ny + 0, 0);
        TVector3D const v0 = fData[i00] + dx * (fData[i10] - fData[i00]) / fXStep;

        size_t const i01 = GetIndex(nx + 0, ny + 1, 0);
        size_t const i11 = GetIndex(nx + 1, ny + 1, 0);
        TVector3D const v1 = fData[i01] + dx * (fData[i11] - fData[i01]) / fXStep;

        return v0 + dy * (v1 - v0) / fYStep;
      }
    case kDIMX_XZ:
      {
        size_t const i00 = GetIndex(nx + 0, 0, nz + 0);
        size_t const i10 = GetIndex(nx + 1, 0, nz + 0);
        TVector3D const v0 = fData[i00] + dx * (fData[i10] - fData[i00]) / fXStep;

        size_t const i01 = GetIndex(nx + 0, 0, nz + 1);
        size_t const i11 = GetIndex(nx + 1, 0, nz + 1);
        TVector3D const v1 = fData[i01] + dx * (fData[i11] - fData[i01]) / fXStep;

        return v0 + dz * (v1 - v0) / fZStep;
      }
    case kDIMX_YZ:
      {
        size_t const i00 = GetIndex(0, ny + 0, nz + 0);
        size_t const i10 = GetIndex(0, ny + 1, nz + 0);
        TVector3D const v0 = fData[i00] + dy * (fData[i10] - fData[i00]) / fYStep;

        size_t const i01 = GetIndex(0, ny + 0, nz + 1);
        size_t const i11 = GetIndex(0, ny + 1, nz + 1);
        TVector3D const v1 = fData[i01] + dy * (fData[i11] - fData[i01]) / fYStep;

        return v0 + dz * (v1 - v0) / fZStep;
      }
      break;
    default:
      throw;
  }

  throw;
}










double TField3DGrid::GetHeaderValue (std::string const& L) const
{
  // Get the value after comment character

  // Easy reading
  std::istringstream S;
  S.str(L);

  // For the two values of interest
  double Value;

  S >> Value;

  // Check read state of input
  if (S.bad()) {
    std::cerr << "ERROR: S is bad" << std::endl;
    throw;
  }

  return Value;
}







double TField3DGrid::GetHeaderValueSRW (std::string const& L, const char CommentChar) const
{
  // Get the value after comment character

  // Easy reading
  std::istringstream S;
  S.str(L);

  // For the two values of interest
  char CC;
  double Value;

  S.get(CC);

  if (CC != CommentChar) {
    std::cerr << "ERROR: bad format in header" << std::endl;
    throw;
  }

  S >> Value;

  // Check read state of input
  if (S.bad()) {
    std::cerr << "ERROR: S is bad" << std::endl;
    throw;
  }

  return Value;
}







void TField3DGrid::ReadFile (std::string const& InFileName, TVector3D const& Rotations, TVector3D const& Translation, char const CommentChar)
{
  // Read file with the best format in the entire world, OSCARSv1.0

  // Open the input file and throw exception if not open
  std::ifstream fi(InFileName);
  if (!fi) {
    std::cerr << "ERROR: cannot open file" << std::endl;
    throw;
  }

  // For reading lines of file
  std::istringstream S;
  std::string L;

  // Initial line for comment
  std::getline(fi, L);


  // Initial X
  std::getline(fi, L);
  double const XStart = GetHeaderValue(L);

  // Step X
  std::getline(fi, L);
  double const XStep = GetHeaderValue(L);

  // Number of points X
  std::getline(fi, L);
  int const NX = (int) GetHeaderValue(L);


  // Initial Y
  std::getline(fi, L);
  double const YStart = GetHeaderValue(L);

  // Step Y
  std::getline(fi, L);
  double const YStep = GetHeaderValue(L);

  // Number of points Y
  std::getline(fi, L);
  int const NY = (int) GetHeaderValue(L);


  // Initial Z
  std::getline(fi, L);
  double const ZStart = GetHeaderValue(L);

  // Step Z
  std::getline(fi, L);
  double const ZStep = GetHeaderValue(L);

  // Number of points Z
  std::getline(fi, L);
  int const NZ = (int) GetHeaderValue(L);


  // Check Number of points is > 0 for all
  if (NX < 1 || NY < 1 || NY < 1) {
    std::cerr << "ERROR: invalid npoints" << std::endl;
    throw;
  }

  // Save position data to object variables
  fNX = NX;
  fNY = NY;
  fNZ = NZ;
  fXStart = XStart;
  fYStart = YStart;
  fZStart = ZStart;
  fXStep  = XStep;
  fYStep  = YStep;
  fZStep  = ZStep;
  fXStop  = fXStart + (fNX - 1) * fXStep;
  fYStop  = fYStart + (fNY - 1) * fYStep;
  fZStop  = fZStart + (fNZ - 1) * fZStep;

  fHasX = NX > 1 ? true : false;
  fHasY = NY > 1 ? true : false;
  fHasZ = NZ > 1 ? true : false;

  if (fHasX && fHasY && fHasZ) {
    fDIMX = kDIMX_XYZ;
  } else if (fHasX && fHasY) {
    fDIMX = kDIMX_XY;
  } else if (fHasX && fHasZ) {
    fDIMX = kDIMX_XZ;
  } else if (fHasY && fHasZ) {
    fDIMX = kDIMX_YZ;
  } else if (fHasX) {
    fDIMX = kDIMX_X;
  } else if (fHasY) {
    fDIMX = kDIMX_Y;
  } else if (fHasZ) {
    fDIMX = kDIMX_Z;
  } else {
    std::cerr << "ERROR: error in file header format" << std::endl;
    throw;
  }

  fXDIM = 0;
  if (fHasX) {
    ++fXDIM;
  }
  if (fHasY) {
    ++fXDIM;
  }
  if (fHasZ) {
    ++fXDIM;
  }

  // Reserve correct number of points in vector (slightly faster)
  fData.reserve(fNX * fNY * fNZ);

  // Temp variables for field
  double fx;
  double fy;
  double fz;


  // Loop over all points
  for (int ix = 0; ix != NX; ++ix) {
    for (int iy = 0; iy != NY; ++iy) {
      for (int iz = 0; iz != NZ; ++iz) {

        // Grab a line from input file
        std::getline(fi, L);

        // Check we did not hit an EOF
        if (fi.eof()) {
          std::cerr << "ERROR: bad input file" << std::endl;
          throw;
        }

        // Read data
        S.clear();
        S.str(L);
        S >> fx >> fy >> fz;

        // Check the stream did not hit an EOF
        if (S.fail()) {
          std::cerr << "ERRROR: input stream bad" << std::endl;
          throw;
        }

        // Push data to storage
        TVector3D F(fx, fy, fz);
        F.RotateSelfXYZ(Rotations);
        fData.push_back(F);
      }
    }
  }

  // Close file
  fi.close();

  // Store Rotations and Translation
  fRotated = Rotations;
  fTranslation = Translation;

  return;
}























void TField3DGrid::ReadFile_SRW (std::string const& InFileName, TVector3D const& Rotations, TVector3D const& Translation, char const CommentChar)
{
  // Read file with SRW field input format

  // Open the input file and throw exception if not open
  std::ifstream fi(InFileName);
  if (!fi) {
    std::cerr << "ERROR: cannot open file" << std::endl;
    throw;
  }

  // For reading lines of file
  std::istringstream S;
  std::string L;

  // Initial line for comment
  std::getline(fi, L);


  // Initial X
  std::getline(fi, L);
  double const XStart = GetHeaderValueSRW(L, CommentChar);

  // Step X
  std::getline(fi, L);
  double const XStep = GetHeaderValueSRW(L, CommentChar);

  // Number of points X
  std::getline(fi, L);
  int const NX = (int) GetHeaderValueSRW(L, CommentChar);


  // Initial Y
  std::getline(fi, L);
  double const YStart = GetHeaderValueSRW(L, CommentChar);

  // Step Y
  std::getline(fi, L);
  double const YStep = GetHeaderValueSRW(L, CommentChar);

  // Number of points Y
  std::getline(fi, L);
  int const NY = (int) GetHeaderValueSRW(L, CommentChar);


  // Initial Z
  std::getline(fi, L);
  double const ZStart = GetHeaderValueSRW(L, CommentChar);

  // Step Z
  std::getline(fi, L);
  double const ZStep = GetHeaderValueSRW(L, CommentChar);

  // Number of points Z
  std::getline(fi, L);
  int const NZ = (int) GetHeaderValueSRW(L, CommentChar);


  // Check Number of points is > 0 for all
  if (NX < 1 || NY < 1 || NY < 1) {
    std::cerr << "ERROR: invalid npoints" << std::endl;
    throw;
  }

  // Save position data to object variables
  fNX = NX;
  fNY = NY;
  fNZ = NZ;
  fXStart = XStart;
  fYStart = YStart;
  fZStart = ZStart;
  fXStep  = XStep;
  fYStep  = YStep;
  fZStep  = ZStep;
  fXStop  = fXStart + (fNX - 1) * fXStep;
  fYStop  = fYStart + (fNY - 1) * fYStep;
  fZStop  = fZStart + (fNZ - 1) * fZStep;

  fHasX = NX > 1 ? true : false;
  fHasY = NY > 1 ? true : false;
  fHasZ = NZ > 1 ? true : false;

  if (fHasX && fHasY && fHasZ) {
    fDIMX = kDIMX_XYZ;
  } else if (fHasX && fHasY) {
    fDIMX = kDIMX_XY;
  } else if (fHasX && fHasZ) {
    fDIMX = kDIMX_XZ;
  } else if (fHasY && fHasZ) {
    fDIMX = kDIMX_YZ;
  } else if (fHasX) {
    fDIMX = kDIMX_X;
  } else if (fHasY) {
    fDIMX = kDIMX_Y;
  } else if (fHasZ) {
    fDIMX = kDIMX_Z;
  } else {
    std::cerr << "ERROR: error in file header format" << std::endl;
    throw;
  }

  fXDIM = 0;
  if (fHasX) {
    ++fXDIM;
  }
  if (fHasY) {
    ++fXDIM;
  }
  if (fHasZ) {
    ++fXDIM;
  }

  // Reserve correct number of points in vector (slightly faster)
  fData.reserve(fNX * fNY * fNZ);

  // Temp variables for field
  double fx;
  double fy;
  double fz;


  // Loop over all points
  for (int ix = 0; ix != NX; ++ix) {
    for (int iy = 0; iy != NY; ++iy) {
      for (int iz = 0; iz != NZ; ++iz) {

        // Grab a line from input file
        std::getline(fi, L);

        // Check we did not hit an EOF
        if (fi.eof()) {
          std::cerr << "ERROR: bad input file" << std::endl;
          throw;
        }

        // Read data
        S.clear();
        S.str(L);
        S >> fx >> fy >> fz;

        // Check the stream did not hit an EOF
        if (S.fail()) {
          std::cerr << "ERRROR: input stream bad" << std::endl;
          throw;
        }

        // Push data to storage
        TVector3D F(fx, fy, fz);
        F.RotateSelfXYZ(Rotations);
        fData.push_back(F);
      }
    }
  }

  // Close file
  fi.close();

  // Store Rotations and Translation
  fRotated = Rotations;
  fTranslation = Translation;

  return;
}








void TField3DGrid::ReadFile_SPECTRA (std::string const& InFileName, TVector3D const& Rotations, TVector3D const& Translation, char const CommentChar)
{
  // Read file with SPECTRA field input format

  // Open the input file and throw exception if not open
  std::ifstream fi(InFileName);
  if (!fi) {
    std::cerr << "ERROR: cannot open file" << std::endl;
    throw;
  }

  // For reading lines of file
  std::istringstream S;
  std::string L;

  // Initial line for comment
  std::getline(fi, L);
  S.str(L);

  // Grab parameters and correct for [mm] -> [m] conversion.
  S >> fXStep >> fYStep >> fZStep >> fNX >> fNY >> fNZ;
  fXStep /= 1000.;
  fYStep /= 1000.;
  fZStep /= 1000.;

  if (S.bad()) {
    throw;
  }


  // Check Number of points is > 0 for all
  if (fNX < 1 || fNY < 1 || fNY < 1) {
    std::cerr << "ERROR: invalid npoints" << std::endl;
    throw;
  }

  // Save position data to object variables
  fXStart = -(fXStep * fNX) / 2;
  fYStart = -(fYStep * fNX) / 2;
  fZStart = -(fZStep * fNX) / 2;
  fXStop  = fXStart + (fNX - 1) * fXStep;
  fYStop  = fYStart + (fNY - 1) * fYStep;
  fZStop  = fZStart + (fNZ - 1) * fZStep;

  fHasX = fNX > 1 ? true : false;
  fHasY = fNY > 1 ? true : false;
  fHasZ = fNZ > 1 ? true : false;

  if (fHasX && fHasY && fHasZ) {
    fDIMX = kDIMX_XYZ;
  } else if (fHasX && fHasY) {
    fDIMX = kDIMX_XY;
  } else if (fHasX && fHasZ) {
    fDIMX = kDIMX_XZ;
  } else if (fHasY && fHasZ) {
    fDIMX = kDIMX_YZ;
  } else if (fHasX) {
    fDIMX = kDIMX_X;
  } else if (fHasY) {
    fDIMX = kDIMX_Y;
  } else if (fHasZ) {
    fDIMX = kDIMX_Z;
  } else {
    std::cerr << "ERROR: error in file header format" << std::endl;
    throw;
  }

  fXDIM = 0;
  if (fHasX) {
    ++fXDIM;
  }
  if (fHasY) {
    ++fXDIM;
  }
  if (fHasZ) {
    ++fXDIM;
  }

  // Reserve correct number of points in vector (slightly faster)
  fData.reserve(fNX * fNY * fNZ);

  // Temp variables for field
  double fx;
  double fy;
  double fz;

  // Loop over all points
  for (int ix = 0; ix != fNX; ++ix) {
    for (int iy = 0; iy != fNY; ++iy) {
      for (int iz = 0; iz != fNZ; ++iz) {

        // Grab a line from input file
        std::getline(fi, L);

        // Check we did not hit an EOF
        if (fi.eof()) {
          std::cerr << "ERROR: bad input file" << std::endl;
          throw;
        }

        // Read data
        S.clear();
        S.str(L);
        S >> fx >> fy >> fz;

        // Check the stream did not hit an EOF
        if (S.fail()) {
          std::cerr << "ERRROR: input stream bad" << std::endl;
          throw;
        }

        // Push data to storage
        TVector3D F(fx, fy, fz);
        F.RotateSelfXYZ(Rotations);
        fData.push_back(F);
      }
    }
  }

  // Close file
  fi.close();

  // Store Rotations and Translation
  fRotated = Rotations;
  fTranslation = Translation;

  return;
}
