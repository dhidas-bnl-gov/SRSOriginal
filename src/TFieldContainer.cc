////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Sep 22 08:19:53 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TFieldContainer.h"




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
