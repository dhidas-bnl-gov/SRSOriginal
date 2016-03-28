////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Mar 25 17:37:39 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <cmath>

#include "TParticleTrajectory2.h"

int TestTParticleTrajectory2 ()
{
  TParticleTrajectory2 P(0.0, 0.001);
  for (int i = 0; i <= 1000; ++i) {
    float ii = i * 0.001;
    P.AddPoint(sin(1 * 2 * 3.14 * ii), ii, ii, ii, ii, ii, ii, ii, ii);
  }

  for (int i = 0; i != 10000; ++i) {
    float ii = i * 0.0001;
    P.GetXAtTime(ii).GetX();
    std::cout << ii << " " << P.GetXAtTime(ii).GetY() << std::endl;
  }


  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TestTParticleTrajectory2();

  return 0;
}
