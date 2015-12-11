////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Dec 10 17:10:50 EST 2015
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TBField1DZ.h"

int TestTBField1DZ (std::string const InFileName)
{
  TBField1DZ TBF;
  TBF.ReadFile(InFileName);
  TBF.Add(1.3, 1.321);
  TBF.Sort();

  TBF.SaveAs("Saved.dat");

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InFileName" << std::endl;
    return 1;
  }

  TestTBField1DZ(argv[1]);

  return 0;
}
