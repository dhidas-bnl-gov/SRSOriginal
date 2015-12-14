#include "TBField1DZ.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

TBField1DZ::TBField1DZ ()
{
  // Default constructor
}



TBField1DZ::TBField1DZ (std::string const& InFileName)
{
  // Constructor.  Reads input file
  this->ReadFile(InFileName);

  // Sort array
  this->Sort();
}



TBField1DZ::~TBField1DZ ()
{
  // Destructor!!
}



bool TBField1DZ::Add (double const Z, double const By)
{
  // Add an element to the array and mark that it is now not necessairly sorted
  fBField.push_back( std::array<double, 2>{ {Z, By} } );
  fIsSorted = false;

  return true;
}



bool TBField1DZ::Sort()
{
  // Sort the field vector by Z, and change IsSorted flag to true
  std::sort(fBField.begin(), fBField.end(), this->CompareBField1DZ);
  fIsSorted = true;

  return true;
}



bool TBField1DZ::IsSorted ()
{
  // Return the state of the IsSorted flag
  return fIsSorted;
}



double TBField1DZ::GetFirstZ ()
{
  // Return the first Z value in the fBField vector.  Make sure it is sorted first
  if (!this->IsSorted()) {
    throw;
  }
  return (*(fBField.begin()))[0];
}



double TBField1DZ::GetLastZ ()
{
  // Return the last Z value in the fBField vector.  Make sure it is sorted first
  if (!this->IsSorted()) {
    throw;
  }
  return (*(fBField.end()))[0];
}




bool TBField1DZ::ReadFile (std::string const& InFileName)
{
  // Read an input file in the format of:
  //   Z By
  // where a line beginning with # is a comment and blank lines are skipped

  // Open the input file.  If !f return a false
  std::ifstream f(InFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: TBField1DZ::ReadFile cannot open file: " << InFileName << std::endl;
    return false;
  }

  // Because we are reading a file this is not necessairly sorted data
  fIsSorted = false;

  // Stream the line for input.  Variables to be filled in look over lines
  std::istringstream Iine;
  double Z, By;

  // Loop over all lines in input file.  Skip if the line is blank or begins with # (a comment)
  for (std::string Line; std::getline(f, Line); ) {
    if (Line == "") {
      continue;
    }

    // Set the streaming line, read it into the doubles
    Iine.clear();
    Iine.str(Line);
    Iine >> Z >> By;

    // Save in field vector
    fBField.push_back(std::array<double, 2>{ {Z, By} });

  }

  // Sort the field from -Z to +Z
  this->Sort();

  return true;
}



bool TBField1DZ::SaveAs (std::string const& OutFileName, std::string const& Comment)
{
  // Save the magnetic field data as a text file with the name given.
  // A comment can be placed at the top in addition to the format comment

  // Print message about saving file
  std::cout << "TBField1DZ::SaveAs saving file as: " << OutFileName << std::endl;

  // Open the output file
  std::ofstream f(OutFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: TBField1DZ::SaveAs cannot open file for writing: " << OutFileName << std::endl;
    return false;
  }

  // Write comment if any, then format comment
  if (Comment != "") {
    f << "# " << Comment << std::endl;
  }
  f << "# [Z in meters] [By in Tesla]" << std::endl;

  // Write BField
  for (std::vector< std::array<double, 2> >::iterator it = fBField.begin(); it != fBField.end(); ++it) {
    f << (*it)[0] << " " << (*it)[1] << "\n";
  }

  // Close file
  f.close();

  // If nothing failed return true
  return true;
}


bool TBField1DZ::CompareBField1DZ (std::array<double, 2> const& A, std::array<double, 2> const& B)
{
  // This function is used for sorting the field in 'Z' cood.  It is a comparison function
  return A[0] < B[0];
}

