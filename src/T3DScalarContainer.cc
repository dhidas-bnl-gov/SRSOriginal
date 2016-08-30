#include "T3DScalarContainer.h"



T3DScalarContainer::T3DScalarContainer ()
{
  // Default constructor
}




T3DScalarContainer::~T3DScalarContainer ()
{
  // Destruction!!!
}




void T3DScalarContainer::AddPoint (TVector3D const& X, double const V)
{
  fValues.push_back( T3DScalar(X, V) );
  fCompensation.push_back(0);
  return;
}




void T3DScalarContainer::AddToPoint (size_t const i, double const V)
{
  // Compensated sum for adding to points

  // Check that the point is within range
  if (i >= fValues.size()) {
    throw std::length_error("T3DScalarContainer::AddtoPoint index out of range");
  }

  double Sum = fValues[i].GetV();
  double y = V - fCompensation[i];
  double t = Sum + y;
  fCompensation[i] = (t - Sum) - y;
  fValues[i].SetV(t);

  return;
}



void T3DScalarContainer::WriteToFileText (std::string const& OutFileName, int const Dimension)
{
  // Write to file in text format
  // If writing to a file, open it and set to scientific output
  std::ofstream of(OutFileName.c_str());
  if (!of.is_open()) {
    throw std::ofstream::failure("cannot open output file");
  }
  of << std::scientific;

  for (size_t i = 0; i != this->GetNPoints(); ++i) {
    TVector3D const& Obs = this->GetPoint(i).GetX();

    if (Dimension == 2) {
      of << Obs.GetX() << " " << Obs.GetY() << " " << this->GetPoint(i).GetV() << "\n";
    } else if (Dimension == 3) {
      of << Obs.GetX() << " " << Obs.GetY() << " " << Obs.GetZ() << " " << this->GetPoint(i).GetV() << "\n";
    } else {
      throw std::out_of_range("incorrect dimensions");
    }

  }

  of.close();

  return;
}




size_t T3DScalarContainer::GetNPoints () const
{
  return fValues.size();
}




T3DScalar const& T3DScalarContainer::GetPoint (size_t const i) const
{
  if (i >= fValues.size()) {
    throw;
  }

  return fValues[i];
}







