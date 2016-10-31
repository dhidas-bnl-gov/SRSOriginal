////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Sep 22 08:19:53 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TFieldContainer.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>



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
  std::cout << "got format " << OutFormat << std::endl;

  // Field for writing
  TVector3D B;

  // Position
  TVector3D X;

  // Open file for output
  std::ofstream of(OutFileName.c_str());
  if (!of.is_open()) {
    throw;
  }

  std::string CommentNoCRLF = Comment;
  std::replace(CommentNoCRLF.begin(), CommentNoCRLF.end(), '\n', ' ');
  std::replace(CommentNoCRLF.begin(), CommentNoCRLF.end(), '\r', ' ');

  // OSCARS format is text by default
  if (OutFormat == "OSCARS") {

    if (CommentNoCRLF == "") {
      of << "# OSCARS format file" << std::endl;
    } else {
      of << "# " << CommentNoCRLF << std::endl;
    }

    double const XStep = (XLim[1] - XLim[0]) / (NX - 1);
    double const YStep = (YLim[1] - YLim[0]) / (NY - 1);
    double const ZStep = (ZLim[1] - ZLim[0]) / (NZ - 1);

    int const Width = 15;
    of << std::setw(Width) << std::left << XLim.GetX() << "  # X Start position"      << std::endl;
    of << std::setw(Width) << std::left << XStep       << "  # Step size in X"        << std::endl;
    of << std::setw(Width) << std::left << NX          << "  # Number of points in X" << std::endl;
    of << std::setw(Width) << std::left << YLim.GetX() << "  # Y Start position"      << std::endl;
    of << std::setw(Width) << std::left << YStep       << "  # Step size in Y"        << std::endl;
    of << std::setw(Width) << std::left << NY          << "  # Number of points in Y" << std::endl;
    of << std::setw(Width) << std::left << ZLim.GetX() << "  # Z Start position"      << std::endl;
    of << std::setw(Width) << std::left << ZStep       << "  # Step size in Z"        << std::endl;
    of << std::setw(Width) << std::left << NZ          << "  # Number of points in Z" << std::endl;

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
  } else if (std::string(OutFormat.begin(), OutFormat.begin() + 8) == "OSCARS1D") {
    std::cout << "Here on export" << std::endl;
  } else if (OutFormat == "SRW") {

    if (CommentNoCRLF == "") {
      of << "# SRW format file generated by OSCARS" << std::endl;
    } else {
      of << "# " << CommentNoCRLF << std::endl;
    }

    double const XStep = (XLim[1] - XLim[0]) / (NX - 1);
    double const YStep = (YLim[1] - YLim[0]) / (NY - 1);
    double const ZStep = (ZLim[1] - ZLim[0]) / (NZ - 1);

    of << "#\t" << XLim.GetX() << "\t# X Start position" << std::endl;
    of << "#\t" << XStep << std::endl;
    of << "#\t" << NX << std::endl;
    of << "#\t" << YLim.GetX() << "  #\tY Start position" << std::endl;
    of << "#\t" << YStep << std::endl;
    of << "#\t" << NY << std::endl;
    of << "#\t" << ZLim.GetX() << "  #\tZ Start position" << std::endl;
    of << "#\t" << ZStep << std::endl;
    of << "#\t" << NZ << std::endl;

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
          of << B.GetX() << "\t" << B.GetY() << "\t" << B.GetZ() << std::endl;
        }
      }
    }
  } else if (OutFormat == "SPECTRA") {

    if (CommentNoCRLF == "") {
      of << "# SPECTRA format file generated by OSCARS" << std::endl;
    } else {
      of << "# " << CommentNoCRLF << std::endl;
    }

    double const XStep = (XLim[1] - XLim[0]) / (NX - 1);
    double const YStep = (YLim[1] - YLim[0]) / (NY - 1);
    double const ZStep = (ZLim[1] - ZLim[0]) / (NZ - 1);

    // Header information
    of << XStep*1000. << " " << YStep*1000. << " " << ZStep*1000. << " " << NX << " " << NY << " " << NZ << std::endl;

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
