////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Mar 30 19:00:32 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include "TBFieldIdeal1D.h"


int del ()
{
  TBFieldIdeal1D B(0.200, 9, 0, 1.0);


  for (double z = -1.10; z <= 1.10; z += 0.001) {
    std::cout << z << " " << B.GetBy(0, 0, z) << std::endl;
  }



  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  del();

  return 0;
}
