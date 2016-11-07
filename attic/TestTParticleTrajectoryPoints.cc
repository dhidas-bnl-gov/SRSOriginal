////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu May 19 13:23:42 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TParticleTrajectoryPoints.h"


int TestTParticleTrajectoryPoints ()
{


  TParticleTrajectoryPoints P(0.001);


  P.AddPoint(1, 2, 3, 4, 5, 6, 7, 8, 9);

  for (int i = 0; i != P.GetNPoints(); ++i) {
    std::cout << "X      " << P.GetX(i)      << std::endl;
    std::cout << "B      " << P.GetB(i)      << std::endl;
    std::cout << "V      " << P.GetV(i)      << std::endl;
    std::cout << "AoverC " << P.GetAoverC(i) << std::endl;
    std::cout << "A      " << P.GetA(i)      << std::endl;
  }

  TVector3D const A(3, 2, 1);

  std::cout << P.GetX(0) + A << std::endl;
  std::cout << A + P.GetX(0) << std::endl;
  std::cout << P.GetV(0) + A << std::endl;
  std::cout << A + P.GetV(0) / 2 << std::endl;

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TestTParticleTrajectoryPoints();

  return 0;
}
