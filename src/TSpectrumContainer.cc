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

  fSpectrumPoints = new std::vector<std::pair<double, double> >();
}




TSpectrumContainer::TSpectrumContainer (std::vector<double> const& v)
{
  // I guess you'll build it yourself by using AddPoint()

  fSpectrumPoints = new std::vector<std::pair<double, double> >();
  fSpectrumPoints->reserve(v.size());

  for (size_t i = 0; i != v.size(); ++i) {
    fSpectrumPoints->push_back( std::make_pair(v[i], 0.0) );
  }

}


TSpectrumContainer::TSpectrumContainer (size_t const N, double const EFirst, double const ELast)
{
  // If you call this with N==1 it will only use EFirst

  fSpectrumPoints = new std::vector<std::pair<double, double> >(N, std::make_pair(0.0, 0.0));

  if (N < 1) {
    throw;
  }

  if (N == 1) {
    (*fSpectrumPoints)[0].first = EFirst;
  }

  for (size_t i = 0; i != fSpectrumPoints->size(); ++i) {
    (*fSpectrumPoints)[i].first = (ELast - EFirst) / (N - 1) * (double) (i + 1);
  }

}




TSpectrumContainer::~TSpectrumContainer ()
{
  delete fSpectrumPoints;
}



void TSpectrumContainer::SetFlux (size_t const i, double const Flux)
{
  if (i < fSpectrumPoints->size()) {
    (*fSpectrumPoints)[i].second = Flux;
    return;
  }

  throw;

  return;
}





void TSpectrumContainer::SetPoint (size_t const i, double const Energy, double const Flux)
{
  (*fSpectrumPoints)[i].first  = Energy;
  (*fSpectrumPoints)[i].second = Flux;

  return;
}




void TSpectrumContainer::AddPoint (double const Energy)
{
  fSpectrumPoints->push_back( std::make_pair(Energy, 0.0) );

  return;
}






double TSpectrumContainer::GetFlux (size_t const i) const
{
  return (*fSpectrumPoints)[i].second;
}






double TSpectrumContainer::GetEnergy (size_t const i) const
{
  return (*fSpectrumPoints)[i].first;
}






double TSpectrumContainer::GetAngularFrequency (size_t const i) const
{
  return (*fSpectrumPoints)[i].first * TSRS::TwoPi() / 4.1357e-15;
}





size_t TSpectrumContainer::GetNPoints () const
{
  return fSpectrumPoints->size();
}




void TSpectrumContainer::SaveToFile (std::string const fn, std::string const Header) const
{
  std::ofstream f(fn.c_str());

  if (Header != "") {
    f << Header << std::endl;
  }


  f << std::scientific;

  for (std::vector<std::pair<double, double> >::iterator it = fSpectrumPoints->begin(); it != fSpectrumPoints->end(); ++it) {
    f << it->first << " " << it->second << std::endl;
  }

  f.close();

  return;
}
