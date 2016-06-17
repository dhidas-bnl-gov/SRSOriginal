////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jun 14 08:16:09 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TBFieldContainer.h"
#include "TBField3DZRegularized.h"

int TestBFC ()
{
  TBFieldContainer B;

  B.AddField( (TBField*)  new TBField3DZRegularized("epu_linear.dat"));


  std::cout << B.GetBy(0, 0, 0) << std::endl;


  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TestBFC();

  return 0;
}
