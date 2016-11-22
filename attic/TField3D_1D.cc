#include "TField3D_1D.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>

TField3D_1D::TField3D_1D ()
{
  // Default constructor
}



TField3D_1D::TField3D_1D (std::string const& InFileName, TVector3D const& Rotation, TVector3D const& Translation, std::vector<double> const& Scaling)
{
  // Constructor.  Reads input file
  this->ReadFile(InFileName, Rotation, Translation, Scaling);

  // Sort array
  this->Sort();
}



TField3D_1D::~TField3D_1D ()
{
  // Destructor!!
}



bool TField3D_1D::Add (double const Z, double const Bx, double const By, double const Bz)
{
  // Add an element to the array and mark that it is now not necessairly sorted
  //fField.push_back( std::array<double, 4>{ {Z, Bx, By, Bz} } );
  // For compatibility with nvcc
  std::array<double, 4> a = { {Z, Bx, By, Bz} };
  fField.push_back(a);

  fIsSorted = false;

  return true;
}



bool TField3D_1D::Sort()
{
  // Sort the field vector by Z, and change IsSorted flag to true
  std::sort(fField.begin(), fField.end(), this->CompareField3D_1D);
  fIsSorted = true;

  return true;
}



bool TField3D_1D::IsSorted ()
{
  // Return the state of the IsSorted flag
  return fIsSorted;
}



double TField3D_1D::GetFirstZ ()
{
  // Return the first Z value in the fField vector.  Make sure it is sorted first
  if (!this->IsSorted()) {
    throw;
  }
  return (*(fField.begin()))[0];
}



double TField3D_1D::GetLastZ ()
{
  // Return the last Z value in the fField vector.  Make sure it is sorted first
  if (!this->IsSorted()) {
    throw;
  }
  return (*(fField.end()))[0];
}




bool TField3D_1D::ReadFile (std::string const& InFileName, TVector3D const& Rotation, TVector3D const& Translation, std::vector<double> const& Scaling)
{
  // Read an input file in the format of:
  //   Z By
  // where a line beginning with # is a comment and blank lines are skipped

  // UPDATE: CRITICAL - implement rotation

  // Open the input file.  If !f, throw
  std::ifstream f(InFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: TField3D_1D::ReadFile cannot open file: " << InFileName << std::endl;
    throw std::ifstream::failure("TField3D_1D::ReadFile cannot open file");
  }

  // Because we are reading a file this is not necessairly sorted data
  fIsSorted = false;

  // Stream the line for input.  Variables to be filled in look over lines
  std::istringstream Iine;
  double Z, Bx, By, Bz;

  // Loop over all lines in input file.  Skip if the line is blank or begins with # (a comment)
  for (std::string Line; std::getline(f, Line); ) {

    // Look for a blank line or comment line and skip if found.  You should never use tab btw.
    size_t FirstChar = Line.find_first_not_of(" \t");
    if (FirstChar == std::string::npos || Line[FirstChar] == '#') {
      continue;
    }

    // Set the streaming line, read it into the doubles
    Iine.clear();
    Iine.str(Line);
    Iine >> Z >> Bx >> By >> Bz;

    // Scaling
    for (size_t is = 0; is != Scaling.size(); ++is) {
      switch (is) {
        case 0:
          Z *= Scaling[0];
          break;
        case 1:
          Bx *= Scaling[1];
          break;
        case 2:
          By *= Scaling[2];
          break;
        case 3:
          Bz *= Scaling[3];
          break;
        default:
          throw std::ifstream::failure("TField3D_1D::ReadFile scaling is incorrect");
      }
    }

    // UPDATE: rotation
    if (Rotation.GetX() != 0 || Rotation.GetY() != 0) {
      throw std::invalid_argument("Rotation in X and Y not allowed here");
    }
    if (Rotation.GetZ() != 0) {
      double const Bx0 = Bx;
      Bx = Bx * cos(Rotation.GetZ()) - By * sin(Rotation.GetZ());
      By = Bx0 * sin(Rotation.GetZ()) + By * cos(Rotation.GetZ());
    }

    // Translation
    Z += Translation.GetZ();


    // Check the stream to see if it is not good
    if (Iine.fail()) {
      std::cerr << "ERROR: TField3D_1D::ReadFile: data format error on this line: " << Line << std::endl;
      throw std::ifstream::failure("TField3D_1D::ReadFile read error");
    }

    // Save in field vector
    //fField.push_back(std::array<double, 4>{ {Z, Bx, By, Bz} });
    // Compatibility with nvcc
    std::array<double, 4> a = { {Z, Bx, By, Bz} };
    fField.push_back(a);

  }

  // Sort the field from -Z to +Z
  this->Sort();

  return true;
}



bool TField3D_1D::SaveAs (std::string const& OutFileName, std::string const& Comment)
{
  // Save the magnetic field data as a text file with the name given.
  // A comment can be placed at the top in addition to the format comment

  // Print message about saving file
  std::cout << "TField3D_1D::SaveAs saving file as: " << OutFileName << std::endl;

  // Open the output file
  std::ofstream f(OutFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: TField3D_1D::SaveAs cannot open file for writing: " << OutFileName << std::endl;
    return false;
  }

  // Write comment if any, then format comment
  if (Comment != "") {
    f << "# " << Comment << std::endl;
  }
  f << "# [Z in meters] [Bx By Bz in Tesla]" << std::endl;

  // Write Field
  for (std::vector< std::array<double, 4> >::iterator it = fField.begin(); it != fField.end(); ++it) {
    f << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << " " << (*it)[3] << "\n";
  }

  // Close file
  f.close();

  // If nothing failed return true
  return true;
}


bool TField3D_1D::Regularize (std::vector<std::array<double, 3> >& oV, double& oFirstZ, double& oLastZ, double& oStepSizeZ, size_t const NPointsPerMeter)
{
  // Regularize the data in fField and save it in oV.  also save the first
  // and last Z values.  The input for resolution is in points per meter.



  // Check to see if we are sorted yet or not and of not, sort
  if (!this->IsSorted()) {
    this->Sort();
  }

  // Check to see there are at least 2 data points
  if (fField.size() < 2) {
    std::cerr << "ERROR: TField3D_1D::Regularize: not enough data points" << std::endl;
    throw std::length_error("TField3D_1D::Regularize not enough points");
  }

  // Grab the first and last Z
  double const First = fField[0][0];
  double const Last  = fField[fField.size() - 1][0];

  // Get the number of points and the step size
  size_t const NPoints  = (Last - First) * NPointsPerMeter;
  double const StepSize = (Last - First) / (double) (NPoints - 1);
  std::cout << "MESSAGE: TField3D_1D::Regularize StepSize: " << StepSize << std::endl;

  // Clear the output vector and allocated space
  oV.clear();
  oV.reserve(NPoints);


  // Variables to hold the slope between two real points and new By (linear interpolated)
  double SlopeX;
  double SlopeY;
  double SlopeZ;
  double NewBx;
  double NewBy;
  double NewBz;

  // I only initialize to avoid compile-time warnings.  These are for finding the
  // adj bins for linear interpolation
  size_t MinBin    = 0;
  size_t AfterBin  = 0;
  size_t BeforeBin = 0;

  // Variable to hold new calculated Z position
  double ThisZ;

  // For each desired point find the bin before and after the desired Z position
  for (size_t i = 0; i != NPoints; ++i) {
    ThisZ = i * StepSize + First;
    for (size_t j = MinBin + 1; j != fField.size(); ++j) {
      if (fField[j][0] > ThisZ) {
        AfterBin  = j;
        BeforeBin = j - 1;
        MinBin = j - 1;
        break;
      }
    }


    // Calculate the By at desired position based on linear interpolation
    SlopeX = (fField[AfterBin][1] - fField[BeforeBin][1]) /  (fField[AfterBin][0] - fField[BeforeBin][0]);
    SlopeY = (fField[AfterBin][2] - fField[BeforeBin][2]) /  (fField[AfterBin][0] - fField[BeforeBin][0]);
    SlopeZ = (fField[AfterBin][3] - fField[BeforeBin][3]) /  (fField[AfterBin][0] - fField[BeforeBin][0]);
    NewBx = fField[BeforeBin][1] + (ThisZ - fField[BeforeBin][0]) * SlopeX;
    NewBy = fField[BeforeBin][2] + (ThisZ - fField[BeforeBin][0]) * SlopeY;
    NewBz = fField[BeforeBin][3] + (ThisZ - fField[BeforeBin][0]) * SlopeZ;

    // Append the new BxByBz to the output vector
    //oV.push_back(std::array<double, 3>{ {NewBx, NewBy, NewBz} });
    std::array<double, 3> a = { {NewBx, NewBy, NewBz} };
    oV.push_back(a);

  }

  // Copy information to output variables
  oFirstZ    = First;
  oLastZ     = Last;
  oStepSizeZ = StepSize;


  return true;
}



bool TField3D_1D::CompareField3D_1D (std::array<double, 4> const& A, std::array<double, 4> const& B)
{
  // This function is used for sorting the field in 'Z' cood.  It is a comparison function
  return A[0] < B[0];
}

