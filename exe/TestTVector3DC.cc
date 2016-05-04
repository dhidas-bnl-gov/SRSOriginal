////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed May  4 14:07:00 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TVector3D.h"
#include "TVector3DC.h"


int TestTVector3DC ()
{

  TVector3D  a(1, 2, 3);
  TVector3DC b(std::complex<double>(1, 1), 0, 0);
  TVector3DC c(std::complex<double>(0, 0), 1, 0);

  std::cout << a << std::endl;
  std::cout << c << std::endl;


  std::cout << b.Cross(c) << std::endl;



  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TestTVector3DC();

  return 0;
}
