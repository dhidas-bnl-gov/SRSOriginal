#include "TBField1DZRegularized.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

TBField1DZRegularized::TBField1DZRegularized ()
{
  // Default constructor

  // Set some defaults
  fZNPointsPerMeter = 0;
}



TBField1DZRegularized::TBField1DZRegularized (std::string const& InFileName)
{
  // Constructor.

  // Set some defaults
  fZNPointsPerMeter = 0;

  // Reads input file, constructs regularized data
  this->ReadFile(InFileName);
}



TBField1DZRegularized::TBField1DZRegularized (std::string const& InFileName, size_t const N)
{
  // Constructor.

  // Set number of points per meter
  fZNPointsPerMeter = N;

  // Reads input file, constructs regularized data
  this->ReadFile(InFileName);
}



TBField1DZRegularized::TBField1DZRegularized (TBField1DZ& BF)
{
  // Constructor.

  // Set number of points per meter
  fZNPointsPerMeter = 0;

  // Regularize
  this->RegularizeField(BF);
}



TBField1DZRegularized::TBField1DZRegularized (TBField1DZ& BF, size_t const N)
{
  // Constructor.

  // Set number of points per meter
  fZNPointsPerMeter = N;

  // Regularize
  this->RegularizeField(BF);
}



TBField1DZRegularized::~TBField1DZRegularized ()
{
  // Destructor!!
}



bool TBField1DZRegularized::ReadFile (std::string const& InFileName)
{
  // This function will take a field in any order, create a TBField1DZ, and
  // store the regularized field while discarding the TBField1DZ.

  TBField1DZ BF(InFileName);
  this->RegularizeField(BF);
  return true;
}



bool TBField1DZRegularized::ReadFileRegularized (std::string const& InFileName)
{
  // Read an input file in the format of:
  //   Z By
  // where a line beginning with # is a comment and blank lines are skipped.
  // Z in this case is ignored since the data is assumed to be regularized.
  // First Z and last Z are used along with NPoints to determine grid.
  // Z should come in increasing order (smallest Z first, largest Z last).
  // It is YOUR responsibility to make sure the input file is already on a 
  // regularized grid for using this function.

  // Open the input file.  If !f return a false
  std::ifstream f(InFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: TBField1DZ::ReadFile cannot open file: " << InFileName << std::endl;
    return false;
  }

  // Stream the line for input.  Variables to be filled in look over lines
  std::istringstream Iine;
  double Z, By;

  // Is this the first real line of input or not?
  bool IsFirst = true;

  // How many points total do we have?  Start with zero
  fZNPoints = 0;

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
    Iine >> Z >> By;

    // Check the stream to see if it is not good
    if (Iine.fail()) {
      std::cerr << "ERROR: TBField1DZRegularized::ReadFileRegularized: data format error on this line: " << Line << std::endl;
      throw;
    }

    // Keep track of the number of points
    ++fZNPoints;

    // This is the first real point
    if (IsFirst) {
      fZFirstPoint = Z;
      IsFirst = false;
    }

    // Save in field vector
    fBField.push_back(By);

  }

  // Set the last Z point
  fZLastPoint = Z;

  // Basic check of inputs
  if (fZNPoints <= 1 || fZLastPoint <= fZFirstPoint) {
    std::cerr << "ERROR: TBField1DZRegularized::ReadFile data format problem in input file" << std::endl;
    throw;
    return false;
  }

  // Set the step size
  fZStepSize = (fZLastPoint - fZFirstPoint) / (fZNPoints - 1.0);

  return true;
  }






bool TBField1DZRegularized::SaveAs (std::string const& OutFileName, std::string const& Comment)
{
  // Save the magnetic field data as a text file with the name given
  // A comment can be placed at the top in addition to the format comment

  // Print message about saving file
  std::cout << "TBField1DZRegularized::SaveAs saving file as: " << OutFileName << std::endl;

  // Open the output file
  std::ofstream f(OutFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: TBField1DZRegularized::SaveAs cannot open file for writing: " << OutFileName << std::endl;
    return false;
  }

  // Write comment if any, then format comment
  if (Comment != "") {
    f << "# " << Comment << std::endl;
  }
  f << "# [Z in meters] [By in Tesla]" << std::endl;

  // Write BField
  double Z;
  for (size_t i = 0; i != fBField.size(); ++i) {
    Z = fZFirstPoint + fZStepSize * (double) i;
    f << Z << " " << fBField[i] << "\n";
  }

  // Close file
  f.close();

  // If nothing failed return true
  return true;
}





double TBField1DZRegularized::GetByAtZ (double const& Z)
{
  // Return the estimated By at a given Z based on grid and linear interpolation

  double const i =  (Z - fZFirstPoint) / fZStepSize;
  int    const j = (int) i;


  //std::cout << "math : " << i << "  " << j << "  " <<(i - (double) j) << "  " << fBField[j+1] << "  " << fBField[j] <<
  //  "  " << (i - (double) j) << "  " << (fBField[j+1] - fBField[j])  << std::endl;

  return fBField[j] + ((fBField[j+1] - fBField[j]) * (i - (double) j));
}



void TBField1DZRegularized::SetZNPointsPerMeter (size_t const N)
{
  // Set the number of points per meter in z direction
  fZNPointsPerMeter = N;

  return;
}




bool TBField1DZRegularized::RegularizeField (TBField1DZ& InField)
{
  if (fZNPointsPerMeter != 0) {
    InField.Regularize(fBField, fZFirstPoint, fZLastPoint, fZStepSize, fZNPointsPerMeter);
  } else {
    InField.Regularize(fBField, fZFirstPoint, fZLastPoint, fZStepSize);
  }
  return true;
}
