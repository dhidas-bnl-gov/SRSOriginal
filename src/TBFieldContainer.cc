////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue May 24 08:59:22 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBFieldContainer.h"




TBFieldContainer::TBFieldContainer ()
{
  // Default constructor
}



TBFieldContainer::TBFieldContainer (TBField* B)
{
  // Construct me with a field
  this->AddField(B);
}



TBFieldContainer::~TBFieldContainer ()
{
  // Default destructor.  I don't own anything I know about
}



void TBFieldContainer::AddField (TBField* B)
{
  // Construct me with a field
  fBFields.push_back(B);
}



double TBFieldContainer::GetBx (double const X, double const Y, double const Z) const
{
  double Sum = 0;

  // Loop over BFields for summing B-fields
  for (std::vector<TBField*>::const_iterator it = fBFields.begin(); it != fBFields.end(); ++it) {
    Sum += (*it)->GetBx(X, Y, Z);
  }

  return Sum;
}



double TBFieldContainer::GetBy (double const X, double const Y, double const Z) const
{
  double Sum = 0;

  // Loop over BFields for summing B-fields
  for (std::vector<TBField*>::const_iterator it = fBFields.begin(); it != fBFields.end(); ++it) {
    Sum += (*it)->GetBy(X, Y, Z);
  }

  return Sum;
}



double TBFieldContainer::GetBz (double const X, double const Y, double const Z) const
{
  double Sum = 0;

  // Loop over BFields for summing B-fields
  for (std::vector<TBField*>::const_iterator it = fBFields.begin(); it != fBFields.end(); ++it) {
    Sum += (*it)->GetBy(X, Y, Z);
  }

  return Sum;
}



TVector3D TBFieldContainer::GetB (double const X, double const Y, double const Z) const
{
  TVector3D Sum(0, 0, 0);

  // Loop over BFields for summing B-fields
  for (std::vector<TBField*>::const_iterator it = fBFields.begin(); it != fBFields.end(); ++it) {
    Sum += (*it)->GetB(X, Y, Z);
  }

  return Sum;
}
