////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu May 19 16:05:08 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include "TSpectrumContainer.h"

#include "TSRS.h"

#include <fstream>



TSpectrumContainer::TSpectrumContainer ()
{
  // I guess you'll build it yourself by using AddPoint()
}




TSpectrumContainer::TSpectrumContainer (std::vector<double> const& V)
{
  // Constructor for an arbitrary list of points
  // V - vector of energy points in [eV]
  this->Init(V);
}


TSpectrumContainer::TSpectrumContainer (size_t const N, double const EFirst, double const ELast)
{
  // Constructor for evenly spaced points in a given energy range.
  // N - Number of points
  // EFirst - First energy point [eV]
  // ELast  - Last energy point [eV]

  this->Init(N, EFirst, ELast);
}




TSpectrumContainer::~TSpectrumContainer ()
{
  // Destruction!
}




void TSpectrumContainer::Init (size_t const N, double const EFirst, double const ELast)
{
  // If you call this with N==1 it will only use EFirst

  // N - Number of points
  // EFirst - First energy point [eV]
  // ELast  - Last energy point [eV]

  // Clear existing data and resize member to the correct size for input
  fSpectrumPoints.clear();
  fSpectrumPoints.resize(N, std::make_pair(0.0, 0.0));


  // If you have zero elements I don't see the point of this
  if (N < 1) {
    throw;
  }

  // If only one point just set it to the 'First' energy
  if (N == 1) {
    fSpectrumPoints[0].first = EFirst;
    return;
  }

  // Set the energy in each element
  for (size_t i = 0; i != fSpectrumPoints.size(); ++i) {
    fSpectrumPoints[i].first = EFirst + (ELast - EFirst) / (N - 1) * (double) (i + 1);
  }

  return;
}




void TSpectrumContainer::Init (std::vector<double> const& V)
{
  // Initialize spectrum with an arbitrary vector of points
  // V - vector of energy points in [eV]

  // Clear existing data and reserve the correct amount for input
  fSpectrumPoints.clear();
  fSpectrumPoints.reserve(V.size());

  // Add each input from V to the internal vector
  for (size_t i = 0; i != V.size(); ++i) {
    fSpectrumPoints.push_back( std::make_pair(V[i], 0.0) );
  }

}



void TSpectrumContainer::SetFlux (size_t const i, double const Flux)
{
  // Set the flux for a given index

  // Simple check
  if (i >= fSpectrumPoints.size()) {
    throw;
  }

  fSpectrumPoints[i].second = Flux;
  return;
}





void TSpectrumContainer::SetPoint (size_t const i, double const Energy, double const Flux)
{
  // Set the energy and flux for a given index

  // I can't decide if I want to be nice and resize or just throw...
  if (i >= fSpectrumPoints.size()) {
    throw;
  }

  fSpectrumPoints[i].first  = Energy;
  fSpectrumPoints[i].second = Flux;

  return;
}




void TSpectrumContainer::AddPoint (double const Energy)
{
  // Add an energy point to the end of the vector.  Default flux value is zero
  fSpectrumPoints.push_back( std::make_pair(Energy, 0.0) );

  return;
}






double TSpectrumContainer::GetFlux (size_t const i) const
{
  // Get flux at a given index
  return fSpectrumPoints[i].second;
}






double TSpectrumContainer::GetEnergy (size_t const i) const
{
  // Get energy at a given index
  return fSpectrumPoints[i].first;
}






double TSpectrumContainer::GetAngularFrequency (size_t const i) const
{
  // Get the angular frequency of this index (from energy)
  // UPDATE: consider removing this function
  return TSRS::EvToAngularFrequency(fSpectrumPoints[i].first);
}





size_t TSpectrumContainer::GetNPoints () const
{
  // Return the number of points in this spectrum
  return fSpectrumPoints.size();
}




void TSpectrumContainer::SaveToFile (std::string const FileName, std::string const Header) const
{
  // Write this spectrum to a file in text format.
  // FileName - File name to write to
  // Header   - Header to print in file

  // Open output file
  std::ofstream f(FileName.c_str());

  // Check if file is open
  // UPDATE: try a more robust check
  if (!f.is_open()) {
    throw;
  }

  // If the header is specified, write one!
  if (Header != "") {
    f << Header << std::endl;
  }

  // Set in scientific mode for printing
  // UPDATE: Could change this to accept c-stype formatting
  f << std::scientific;

  // Loop over spectrum and print to file
  for (std::vector<std::pair<double, double> >::const_iterator it = fSpectrumPoints.begin(); it != fSpectrumPoints.end(); ++it) {
    f << it->first << " " << it->second << std::endl;
  }

  // Close file
  f.close();

  return;
}






void TSpectrumContainer::SaveToFileBinary (std::string const FileName, std::string const Header) const
{
  // UPDATE: Actually make this a binary output
  // Write this spectrum to a file in text format.
  // FileName - File name to write to
  // Header   - Header to print in file

  throw;

  // Open output file
  std::ofstream f(FileName.c_str());

  // Check if file is open
  // UPDATE: try a more robust check
  if (!f.is_open()) {
    throw;
  }

  // If the header is specified, write one!
  if (Header != "") {
    f << Header << std::endl;
  }

  // Set in scientific mode for printing
  // UPDATE: Could change this to accept c-stype formatting
  f << std::scientific;

  // Loop over spectrum and print to file
  for (std::vector<std::pair<double, double> >::const_iterator it = fSpectrumPoints.begin(); it != fSpectrumPoints.end(); ++it) {
    f << it->first << " " << it->second << std::endl;
  }

  // Close file
  f.close();

  return;
}
