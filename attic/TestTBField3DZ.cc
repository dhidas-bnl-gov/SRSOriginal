////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Feb 17 18:44:01 EST 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TBField3DZ.h"

int TestTBField3DZ (std::string const InFileName)
{
  TBField3DZ TBF;
  TBF.ReadFile(InFileName);
  //TBF.Add(1.3, 1.321);
  TBF.Sort();

  TBF.SaveAs("Saved.dat");



  std::vector<std::array<double, 3> > A;
  double FirstZ, LastZ, StepSize;

  TBF.Regularize(A, FirstZ, LastZ, StepSize, 10000);


  for (size_t i = 0; i != A.size(); ++i) {
    std::cout << FirstZ + StepSize * (double) i << "   " << A[i][0] << "   " << A[i][1] << "   " << A[i][2] << "\n";
  }

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InFileName]" << std::endl;
    return 1;
  }

  TestTBField3DZ(argv[1]);

  return 0;
}
