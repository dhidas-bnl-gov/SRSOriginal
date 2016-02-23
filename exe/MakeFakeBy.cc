////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Feb 16 11:55:52 EST 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <string>
#include <fstream>


int MakeFakeBy (std::string const OutFileName)
{
  // Open output file
  std::ofstream f(OutFileName.c_str());
  if (!f.is_open()) {
    std::cerr << "ERROR: Cannot open file: " << OutFileName << std::endl;;
  }


  float const XStart = -0.9;
  float const XStop  =  0.9;

  int const NPoints = 1000;

  float const StepSize = (XStop - XStart) / (NPoints - 1);

  float X, By;
  for (int i = 0; i != NPoints; ++i) {
    X  = XStart + StepSize * i;

    By = 0;

    if (X < 0.0) {
      By = 0;
    } else if (X < 0.1) {
      By = 1;
    } else if (X < 0.3) {
      By = -1;
    } else if (X < 0.4) {
      By = 1;
    } else {
      By = 0;
    }

    f << X << " " << By << std::endl;
  }




  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [OutFileName]" << std::endl;
    return 1;
  }

  MakeFakeBy(argv[1]);

  return 0;
}
