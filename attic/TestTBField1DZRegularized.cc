////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Dec 11 16:53:32 EST 2015
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TBField1DZRegularized.h"

int TestTBField1DZRegularized (std::string const InFileName)
{
  TBField1DZRegularized TBF;
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

  TestTBField1DZRegularized(argv[1]);

  return 0;
}
