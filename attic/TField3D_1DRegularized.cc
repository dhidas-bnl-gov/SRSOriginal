#include "TField3D_1DRegularized.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

TField3D_1DRegularized::TField3D_1DRegularized ()
{
  // Default constructor

  // Set some defaults
  fZNPointsPerMeter = 0;
}



TField3D_1DRegularized::TField3D_1DRegularized (std::string const& InFileName)
{
  // Constructor.

  // Set some defaults
  fZNPointsPerMeter = 0;

  // Reads input file, constructs regularized data
  this->ReadFile(InFileName);
}



TField3D_1DRegularized::TField3D_1DRegularized (std::string const& InFileName, TVector3D const& Rotation, TVector3D const& Translation, std::vector<double> const& Scaling)
{
  // Constructor.

  // Set some defaults
  fZNPointsPerMeter = 0;

  // Reads input file, constructs regularized data
  this->ReadFile(InFileName, Rotation, Translation, Scaling);
}



TField3D_1DRegularized::TField3D_1DRegularized (std::string const& InFileName, size_t const N)
{
  // Constructor.

  // Set number of points per meter
  fZNPointsPerMeter = N;

  // Reads input file, constructs regularized data
  this->ReadFile(InFileName);
}



TField3D_1DRegularized::TField3D_1DRegularized (TField3D_1D& F)
{
  // Constructor.

  // Set number of points per meter
  fZNPointsPerMeter = 0;

  // Regularize
  this->RegularizeField(F);
}



TField3D_1DRegularized::TField3D_1DRegularized (TField3D_1D& F, size_t const N)
{
  // Constructor.

  // Set number of points per meter
  fZNPointsPerMeter = N;

  // Regularize
  this->RegularizeField(F);
}



TField3D_1DRegularized::~TField3D_1DRegularized ()
{
  // Destructor!!
}



bool TField3D_1DRegularized::ReadFile (std::string const& InFileName)
{
  // This function will take a field in any order, create a TField3D_1D, and
  // store the regularized field while discarding the TField3D_1D.

  TField3D_1D F(InFileName);
  this->RegularizeField(F);
  return true;
}



bool TField3D_1DRegularized::ReadFile (std::string const& InFileName, TVector3D const& Rotation, TVector3D const& Translation, std::vector<double> const& Scaling)
{
  // UPDATE: Comment function

  TField3D_1D F(InFileName, Rotation, Translation, Scaling);
  this->RegularizeField(F);
  return true;
}



bool TField3D_1DRegularized::ReadFileRegularized (std::string const& InFileName)
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
    std::cerr << "ERROR: TField3D_1D::ReadFile cannot open file: " << InFileName << std::endl;
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
      std::cerr << "ERROR: TField3D_1DRegularized::ReadFileRegularized: data format error on this line: " << Line << std::endl;
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
    //fField.push_back({ {Bx, By, Bz} });
    // For nvcc compatibility
    std::array<double, 3> myarray = { {Bx, By, Bz} };
    fField.push_back(myarray);


  }

  // Set the last Z point
  fZLastPoint = Z;

  // Basic check of inputs
  if (fZNPoints <= 1 || fZLastPoint <= fZFirstPoint) {
    std::cerr << "ERROR: TField3D_1DRegularized::ReadFile data format problem in input file" << std::endl;
    throw;
    return false;
  }

  // Set the step size
  fZStepSize = (fZLastPoint - fZFirstPoint) / (fZNPoints - 1.0);

  return true;
  }






bool TField3D_1DRegularized::SaveAs (std::string const& OutFileName, std::string const& Comment)
{
  // Save the magnetic field data as a text file with the name given
  // A comment can be placed at the top in addition to the format comment

  // Print message about saving file
  std::cout << "TField3D_1DRegularized::SaveAs saving file as: " << OutFileName << std::endl;

  // Open the output file
  std::ofstream f(OutFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: TField3D_1DRegularized::SaveAs cannot open file for writing: " << OutFileName << std::endl;
    return false;
  }

  // Write comment if any, then format comment
  if (Comment != "") {
    f << "# " << Comment << std::endl;
  }
  f << "# [Z in meters] [By in Tesla]" << std::endl;

  // Write Field
  double Z;
  for (size_t i = 0; i != fField.size(); ++i) {
    Z = fZFirstPoint + fZStepSize * (double) i;
    f << Z << " " << fField[i][0] << " " << fField[i][1]<< " " << fField[i][2] << "\n";
  }

  // Close file
  f.close();

  // If nothing failed return true
  return true;
}




double TField3D_1DRegularized::GetFx (double const X, double const Y, double const Z) const
{
  return this->GetFxAtZ(Z);
}



double TField3D_1DRegularized::GetFy (double const X, double const Y, double const Z) const
{
  return this->GetFyAtZ(Z);
}




double TField3D_1DRegularized::GetFz (double const X, double const Y, double const Z) const
{
  return this->GetFzAtZ(Z);
}




TVector3D TField3D_1DRegularized::GetF (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z));
}




TVector3D TField3D_1DRegularized::GetF (TVector3D const& P) const
{
  return TVector3D(this->GetFxAtZ(P.GetZ()), this->GetFyAtZ(P.GetZ()), this->GetFzAtZ(P.GetZ()));
}








double TField3D_1DRegularized::GetFxAtZ (double const Z) const
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
  return fField[j][0] + ((fField[j+1][0] - fField[j][0]) * (i - (double) j));
}



double TField3D_1DRegularized::GetFyAtZ (double const Z) const
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
  return fField[j][1] + ((fField[j+1][1] - fField[j][1]) * (i - (double) j));
}



double TField3D_1DRegularized::GetFzAtZ (double const Z) const
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
  return fField[j][2] + ((fField[j+1][2] - fField[j][2]) * (i - (double) j));
}



void TField3D_1DRegularized::SetZNPointsPerMeter (size_t const N)
{
  // Set the number of points per meter in z direction
  fZNPointsPerMeter = N;

  return;
}




bool TField3D_1DRegularized::RegularizeField (TField3D_1D& InField)
{
  // If the number of points per meter is specified use it, otherwise use the default in TField3D_1D
  if (fZNPointsPerMeter != 0) {
    InField.Regularize(fField, fZFirstPoint, fZLastPoint, fZStepSize, fZNPointsPerMeter);
  } else {
    InField.Regularize(fField, fZFirstPoint, fZLastPoint, fZStepSize);
  }
  return true;
}
