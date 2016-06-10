////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Jun 10 09:18:27 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TParticleBeamContainer.h"

int TestParticleBeamContainer ()
{

  TParticleBeamContainer PBC;






  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TestParticleBeamContainer();

  return 0;
}
