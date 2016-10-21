////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Sep 22 08:19:53 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TFieldContainer.h"

#include <iostream>
#include <fstream>



TFieldContainer::TFieldContainer ()
{
  // Default constructor
}



TFieldContainer::TFieldContainer (TField* F)
{
  // Construct me with a field
  this->AddField(F);
}



TFieldContainer::~TFieldContainer ()
{
  // Default destructor.  I own everything you have passed me.  Make no mistake there!!
  this->Clear();
}



void TFieldContainer::AddField (TField* F)
{
  // Construct me with a field
  fFields.push_back(F);
}



double TFieldContainer::GetFx (double const X, double const Y, double const Z) const
{
  double Sum = 0;

  // Loop over Fields for summing fields
  for (std::vector<TField*>::const_iterator it = fFields.begin(); it != fFields.end(); ++it) {
    Sum += (*it)->GetFx(X, Y, Z);
  }

  return Sum;
}



double TFieldContainer::GetFy (double const X, double const Y, double const Z) const
{
  double Sum = 0;

  // Loop over Fields for summing fields
  for (std::vector<TField*>::const_iterator it = fFields.begin(); it != fFields.end(); ++it) {
    Sum += (*it)->GetFy(X, Y, Z);
  }

  return Sum;
}



double TFieldContainer::GetFz (double const X, double const Y, double const Z) const
{
  double Sum = 0;

  // Loop over Fields for summing fields
  for (std::vector<TField*>::const_iterator it = fFields.begin(); it != fFields.end(); ++it) {
    Sum += (*it)->GetFz(X, Y, Z);
  }

  return Sum;
}



TVector3D TFieldContainer::GetF (double const X, double const Y, double const Z) const
{
  TVector3D Sum(0, 0, 0);

  // Loop over Fields for summing fields
  for (std::vector<TField*>::const_iterator it = fFields.begin(); it != fFields.end(); ++it) {
    Sum += (*it)->GetF(X, Y, Z);
  }

  return Sum;
}




TVector3D TFieldContainer::GetF (TVector3D const& X) const
{
  TVector3D Sum(0, 0, 0);

  // Loop over Fields for summing fields
  for (std::vector<TField*>::const_iterator it = fFields.begin(); it != fFields.end(); ++it) {
    Sum += (*it)->GetF(X);
  }

  return Sum;
}




size_t TFieldContainer::GetNFields () const
{
  // Return the number of fields input
  return fFields.size();
}




void TFieldContainer::Clear ()
{
  for (std::vector<TField*>::iterator it = fFields.begin(); it != fFields.end(); ++it) {
    if (*it != 0x0) {
      delete *it;
    }
  }

  fFields.clear();

  return;
}





void TFieldContainer::WriteToFile (std::string const& OutFileName, std::string const& OutFormat, TVector2D const& XLim, int const NX, TVector2D const& YLim, int const NY, TVector2D const& ZLim, int const NZ, std::string const Comment)
{
  // Write the magnetic field in a given range to an output file of the chosen format

  // Field for writing
  TVector3D B;

  // Position
  TVector3D X;

  // Open file for output
  std::ofstream of(OutFileName.c_str());
  if (!of.is_open()) {
    throw;
  }

  // OSCARS format is text by default
  if (OutFormat == "OSCARS") {

    // UPDATE: Remove LF and CR from Comment
    if (Comment == "") {
      of << "# OSCARS format file" << std::endl;
    } else {
      of << "# " << Comment << std::endl;
    }

    double const XStep = (XLim[1] - XLim[0]) / (NX - 1);
    double const YStep = (YLim[1] - YLim[0]) / (NY - 1);
    double const ZStep = (ZLim[1] - ZLim[0]) / (NZ - 1);

    of << XLim.GetX() << "  # X Start position" << std::endl;
    of << XStep << std::endl;
    of << NX << std::endl;
    of << YLim.GetX() << "  # Y Start position" << std::endl;
    of << YStep << std::endl;
    of << NY << std::endl;
    of << ZLim.GetX() << "  # Z Start position" << std::endl;
    of << ZStep << std::endl;
    of << NZ << std::endl;

    // Set output format
    of << std::scientific;

    // Loop over all points and output
    for (int i = 0; i < NX; ++i) {
      for (int j = 0; j < NY; ++j) {
        for (int k = 0; k < NZ; ++k) {

          // Set current position
          X.SetXYZ(XLim[0] + XStep * i, YLim[0] + YStep * j, ZLim[0] + ZStep * k);

          // Get B Field
          B = this->GetF(X);

          // Print field to file
          of << B.GetX() << " " << B.GetY() << " " << B.GetZ() << std::endl;
        }
      }
    }
  }


  // Close output file
  of.close();

  return;
}
