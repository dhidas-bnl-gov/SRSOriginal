#include "TBField3DZRegularized.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

TBField3DZRegularized::TBField3DZRegularized ()
{
  // Default constructor

  // Set some defaults
  fZNPointsPerMeter = 0;
}



TBField3DZRegularized::TBField3DZRegularized (std::string const& InFileName)
{
  // Constructor.

  // Set some defaults
  fZNPointsPerMeter = 0;

  // Reads input file, constructs regularized data
  this->ReadFile(InFileName);
}



TBField3DZRegularized::TBField3DZRegularized (std::string const& InFileName, size_t const N)
{
  // Constructor.

  // Set number of points per meter
  fZNPointsPerMeter = N;

  // Reads input file, constructs regularized data
  this->ReadFile(InFileName);
}



TBField3DZRegularized::TBField3DZRegularized (TBField3DZ& BF)
{
  // Constructor.

  // Set number of points per meter
  fZNPointsPerMeter = 0;

  // Regularize
  this->RegularizeField(BF);
}



TBField3DZRegularized::TBField3DZRegularized (TBField3DZ& BF, size_t const N)
{
  // Constructor.

  // Set number of points per meter
  fZNPointsPerMeter = N;

  // Regularize
  this->RegularizeField(BF);
}



TBField3DZRegularized::~TBField3DZRegularized ()
{
  // Destructor!!
}



bool TBField3DZRegularized::ReadFile (std::string const& InFileName)
{
  // This function will take a field in any order, create a TBField3DZ, and
  // store the regularized field while discarding the TBField3DZ.

  TBField3DZ BF(InFileName);
  this->RegularizeField(BF);
  return true;
}



bool TBField3DZRegularized::ReadFileRegularized (std::string const& InFileName)
{
  // Read an input file in the format of:
  //   Z Bx By Bz
  // where a line beginning with # is a comment and blank lines are skipped.
  // Z in this case is ignored since the data is assumed to be regularized.
  // First Z and last Z are used along with NPoints to determine grid.
  // Z should come in increasing order (smallest Z first, largest Z last).
  // It is YOUR responsibility to make sure the input file is already on a 
  // regularized grid for using this function.

  // Open the input file.  If !f return a false
  std::ifstream f(InFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: TBField3DZ::ReadFile cannot open file: " << InFileName << std::endl;
    return false;
  }

  // Stream the line for input.  Variables to be filled in look over lines
  std::istringstream Iine;
  double Z, Bx, By, Bz;

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
    Iine >> Z >> Bx >> By >> Bz;


    // Check the stream to see if it is not good
    if (Iine.fail()) {
      std::cerr << "ERROR: TBField3DZRegularized::ReadFileRegularized: data format error on this line: " << Line << std::endl;
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
    fBField.push_back({ {Bx, By, Bz} });

  }

  // Set the last Z point
  fZLastPoint = Z;

  // Basic check of inputs
  if (fZNPoints <= 1 || fZLastPoint <= fZFirstPoint) {
    std::cerr << "ERROR: TBField3DZRegularized::ReadFile data format problem in input file" << std::endl;
    throw;
    return false;
  }

  // Set the step size
  fZStepSize = (fZLastPoint - fZFirstPoint) / (fZNPoints - 1.0);

  return true;
  }






bool TBField3DZRegularized::SaveAs (std::string const& OutFileName, std::string const& Comment)
{
  // Save the magnetic field data as a text file with the name given
  // A comment can be placed at the top in addition to the format comment

  // Print message about saving file
  std::cout << "TBField3DZRegularized::SaveAs saving file as: " << OutFileName << std::endl;

  // Open the output file
  std::ofstream f(OutFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: TBField3DZRegularized::SaveAs cannot open file for writing: " << OutFileName << std::endl;
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
    f << Z << " " << fBField[i][0] << " " << fBField[i][1]<< " " << fBField[i][2] << "\n";
  }

  // Close file
  f.close();

  // If nothing failed return true
  return true;
}




double TBField3DZRegularized::GetBx (double const X, double const Y, double const Z) const
{
  return this->GetBxAtZ(Z);
}



double TBField3DZRegularized::GetBy (double const X, double const Y, double const Z) const
{
  return this->GetByAtZ(Z);
}




double TBField3DZRegularized::GetBz (double const X, double const Y, double const Z) const
{
  return this->GetBzAtZ(Z);
}




TVector3D TBField3DZRegularized::GetB (double const X, double const Y, double const Z) const
{
  return TVector3D(this->GetBxAtZ(Z), this->GetByAtZ(Z), this->GetBzAtZ(Z));
}








double TBField3DZRegularized::GetBxAtZ (double const Z) const
{
  // Return the estimated Bx at a given Z based on grid and linear interpolation.
  // If the requested Z position is outside of fZFirstPoint and fZLastPoint zero is
  // returned.

  // If requested Z is outside of the data assume zero field
  if (Z < fZFirstPoint or Z > fZLastPoint) {
    return 0;
  }

  // Where in the data array is this value
  double const i =  (Z - fZFirstPoint) / fZStepSize;
  int    const j = (int) i;

  // Return linear estimate of the field based on the two points on either side
  return fBField[j][0] + ((fBField[j+1][0] - fBField[j][0]) * (i - (double) j));
}



double TBField3DZRegularized::GetByAtZ (double const Z) const
{
  // Return the estimated By at a given Z based on grid and linear interpolation.
  // If the requested Z position is outside of fZFirstPoint and fZLastPoint zero is
  // returned.

  // If requested Z is outside of the data assume zero field
  if (Z < fZFirstPoint or Z > fZLastPoint) {
    return 0;
  }

  // Where in the data array is this value
  double const i =  (Z - fZFirstPoint) / fZStepSize;
  int    const j = (int) i;

  // Return linear estimate of the field based on the two points on either side
  return fBField[j][1] + ((fBField[j+1][1] - fBField[j][1]) * (i - (double) j));
}



double TBField3DZRegularized::GetBzAtZ (double const Z) const
{
  // Return the estimated Bz at a given Z based on grid and linear interpolation.
  // If the requested Z position is outside of fZFirstPoint and fZLastPoint zero is
  // returned.

  // If requested Z is outside of the data assume zero field
  if (Z < fZFirstPoint or Z > fZLastPoint) {
    return 0;
  }

  // Where in the data array is this value
  double const i =  (Z - fZFirstPoint) / fZStepSize;
  int    const j = (int) i;

  // Return linear estimate of the field based on the two points on either side
  return fBField[j][2] + ((fBField[j+1][2] - fBField[j][2]) * (i - (double) j));
}



void TBField3DZRegularized::SetZNPointsPerMeter (size_t const N)
{
  // Set the number of points per meter in z direction
  fZNPointsPerMeter = N;

  return;
}




bool TBField3DZRegularized::RegularizeField (TBField3DZ& InField)
{
  // If the number of points per meter is specified use it, otherwise use the default in TBField3DZ
  if (fZNPointsPerMeter != 0) {
    InField.Regularize(fBField, fZFirstPoint, fZLastPoint, fZStepSize, fZNPointsPerMeter);
  } else {
    InField.Regularize(fBField, fZFirstPoint, fZLastPoint, fZStepSize);
  }
  return true;
}
