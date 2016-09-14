////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jul  7 18:57:49 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TSurfacePoints_Rectangle.h"
#include "TSRS.h"

#include <iostream>


int TestRect ()
{
  TSurfacePoints_Rectangle S;
  S.Init("ZX", 2, 2, 1.5, 2.5, TVector3D(0, 0, 0), TVector3D(0, 0, 0), 0);

  std::cout << S.GetNPoints() << std::endl;
  std::cout << S.GetXYZ(0) << std::endl;
  std::cout << S.GetX1(0) << std::endl;
  std::cout << S.GetX2(0) << std::endl;

  for (size_t i = 0; i != S.GetNPoints(); ++i) {
    std::cout << i << "  " << S.GetXYZ(i) << "  " << S.GetPoint(i).GetNormal()<< std::endl;
  }

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TestRect();

  return 0;
}
