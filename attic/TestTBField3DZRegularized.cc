////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Feb 18 17:58:00 EST 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TBField3DZRegularized.h"

int TestTBField3DZRegularized (std::string const InFileName)
{
  TBField3DZRegularized TBF;
  TBF.ReadFileRegularized(InFileName);
  TBF.SaveAs("Saved2.dat");

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InFileName" << std::endl;
    return 1;
  }

  TestTBField3DZRegularized(argv[1]);

  return 0;
}
