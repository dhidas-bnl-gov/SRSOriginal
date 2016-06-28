////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jun 23 08:52:35 EDT 2016
//
// Intended to be an interface to Mathematica for SRS
// including a vector of SRS class objects
//
////////////////////////////////////////////////////////////////////
#include "wstp.h"

#include "SRS.h"
#include "SRS_Mathematica.h"

#include <vector>

std::vector<SRS*> fSRSVector;

int SRS_Init ()
{
  fSRSVector.push_back(new SRS());
  return fSRSVector.size() - 1;
}


int SRS_Delete (int const i)
{
  if (i < fSRSVector.size() && fSRSVector[i] != 0x0) {
    delete fSRSVector[i];
    fSRSVector[i] = 0x0;
    return -1;
  }

  return -1;
}



int SRS_DeleteAll ()
{
  for (size_t i = 0; i != fSRSVector.size(); ++i) {
    if (fSRSVector[i] != 0x0) {
      delete fSRSVector[i];
      fSRSVector[i] = 0x0;
    }
  }

  fSRSVector.clear();

  return -1;
}



void SRS_AddMagneticField (int const i, const char* FileName, const char* Format)
{
  if (i < fSRSVector.size() && fSRSVector[i] != 0x0) {
    fSRSVector[i]->AddMagneticField(FileName, Format);
  }

  return;
}



void SRS_GetB (int const i, double X[], long N)
{
  if (N != 3) {
    throw;
  }


  if (i >= fSRSVector.size() || fSRSVector[i] == 0x0) {
    throw;
  }

  TVector3D B = fSRSVector[i]->GetB(X[0], X[1], X[2]);

  double BReturn[3];
  BReturn[0] = B.GetX();
  BReturn[1] = B.GetY();
  BReturn[2] = B.GetZ();

  WSPutReal64List(stdlink, BReturn, 3);

  return;
}



double SRS_GetBx (int const i, double X[], long N)
{
  if (N != 3) {
    throw;
  }

  if (i >= fSRSVector.size() || fSRSVector[i] == 0x0) {
    throw;
  }

  return fSRSVector[i]->GetBx(X[0], X[1], X[2]);
}



double SRS_GetBy (int const i, double X[], long N)
{
  if (N != 3) {
    throw;
  }

  if (i >= fSRSVector.size() || fSRSVector[i] == 0x0) {
    throw;
  }

  return fSRSVector[i]->GetBy(X[0], X[1], X[2]);
}



double SRS_GetBz (int const i, double X[], long N)
{
  if (N != 3) {
    throw;
  }

  if (i >= fSRSVector.size() || fSRSVector[i] == 0x0) {
    throw;
  }

  return fSRSVector[i]->GetBz(X[0], X[1], X[2]);
}
