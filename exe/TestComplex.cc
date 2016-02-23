////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Feb 17 15:39:27 EST 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <complex>


int TestComplex ()
{

  std::complex<double> c(1, 1);
  std::cout << std::arg(c) << std::endl;
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TestComplex();

  return 0;
}
