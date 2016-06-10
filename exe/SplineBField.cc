////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue May 10 08:56:49 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TROOT.h"

int SplineBField(char* fin, char* fon)
{
  std::ifstream fi(fin);

  double Z, Bx, By, Bz;
  TGraph gx, gy, gz;

  double FirstZ = 99999;
  double LastZ;

  while (true) {
    fi >> Z >> Bx >> By >> Bz;
    if (fi.eof()) {
      break;
    }

    if (Z < FirstZ) {
      FirstZ = Z;
    }
    LastZ = Z;

    gx.Set(gx.GetN()+1);
    gx.SetPoint(gx.GetN()-1, Z, Bx);

    gy.Set(gy.GetN()+1);
    gy.SetPoint(gy.GetN()-1, Z, By);

    gz.Set(gz.GetN()+1);
    gz.SetPoint(gz.GetN()-1, Z, Bz);
  }

  fi.close();

  std::ofstream fo(fon);

 
  std::cout << gx.GetN() << std::endl; 
  std::cout << gy.GetN() << std::endl; 
  std::cout << gz.GetN() << std::endl; 

  int const NewN = 10001;

  double const StepSize = (LastZ - FirstZ) / (double) (NewN - 1);


  TGraph gnx(NewN);
  TGraph gny(NewN);
  TGraph gnz(NewN);

  for (size_t i = 0; i != NewN; ++i) {
    Z = FirstZ + StepSize * (double) i;

    Bx = gx.Eval(Z);
    By = gy.Eval(Z);
    Bz = gz.Eval(Z);

    gnx.SetPoint(i, Z, Bx);
    gny.SetPoint(i, Z, By);
    gnz.SetPoint(i, Z, Bz);
    
    fo << Z << " " << Bx << " " << By << " " << Bz << std::endl;

  }

  gnx.SetMarkerColor(2);
  gny.SetMarkerColor(2);
  gnz.SetMarkerColor(2);

  TCanvas c;
  c.cd();

  gnx.Draw("AP");
  gx.Draw("PLsame");
  c.SaveAs("Compare_Bx.pdf");

  gny.Draw("AP");
  gy.Draw("Psame");
  c.SaveAs("Compare_By.pdf");

  gnz.Draw("AP");
  gz.Draw("Psame");
  c.SaveAs("Compare_Bz.pdf");



  fo.close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [InFile] [OutFile]" << std::endl;
    return 1;
  }

  SplineBField(argv[1], argv[2]);

  return 0;
}
