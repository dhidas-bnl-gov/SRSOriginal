////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Apr  6 15:52:04 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TROOT.h"

int RetimeEField (char* fn)
{
  std::ifstream fi(fn);

  double t, x;
  TGraph g;
  std::vector<double> T, X;
  while (true) {
    fi >> t >> x;
    if (fi.eof()) {
      break;
    }


    g.Set(g.GetN()+1);
    g.SetPoint(g.GetN(), t, x);
    T.push_back(t);
    X.push_back(x);
  }


  

  int const NewN = (int) T.size() * 7;


  double const StepSize = (T[T.size() - 1] - T[0]) / (double) (NewN - 1);


  TCanvas c;
  c.cd();
  g.Draw("AP");

  c.SaveAs("del.png");

  double time, val;
  for (size_t i = 0; i != NewN; ++i) {
    time = T[0] + StepSize * (double) i;
    val = g.Eval(time);

    printf("%40.30E  %40.30E\n", time, val);
  }

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InFile]" << std::endl;
    return 1;
  }

  RetimeEField(argv[1]);

  return 0;
}
