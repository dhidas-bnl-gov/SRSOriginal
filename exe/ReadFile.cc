////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu May 26 17:01:02 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>


bool SortVector1 (std::vector<double> const& L, std::vector<double> const& R)
{
  if (L[0] > R[0]) {
    return false;
  }

  return true;
}
bool SortVector2 (std::vector<double> const& L, std::vector<double> const& R)
{
  if (L[0] > R[0]) {
    return false;
  } else if (L[0] == R[0]) {
    if (L[1] > R[1]) {
      return false;
    } else {
      return true;
    }
  }

  return true;
}
bool SortVector3 (std::vector<double> const& L, std::vector<double> const& R)
{
  if (L[0] > R[0]) {
    return false;
  } else if (L[0] == R[0]) {
    if (L[1] > R[1]) {
      return false;
    } else if (L[1] == R[1]) {
      if (L[2] >= R[2]) {
        return false;
      } else {
        return true;
      }
    }
  } else {
    return true;
  }

  return true;
}







int ReadFile (std::string const InFileName, std::string const Format)
{


  std::ifstream fi(InFileName.c_str());
  if (!fi) {
    std::cerr << "ERROR: cannot read file: " << InFileName << std::endl;
    throw;
  }

  std::vector<bool> HasXB(6, false);

  std::vector<int> Order(6, -1);

  std::istringstream s;
  s.str(Format);

  std::string c;
  int i = 0;
  int XDIM = 0;
  int BDIM = 0;
  while (s >> c) {


    if (c == "X") {
      Order[0] = i;
      ++XDIM;
    } else if (c == "Y") {
      Order[1] = i;
      ++XDIM;
    } else if (c == "Z") {
      Order[2] = i;
      ++XDIM;
    } else if (c == "Bx") {
      Order[3] = i;
      ++BDIM;
    } else if (c == "By") {
      Order[4] = i;
      ++BDIM;
    } else if (c == "Bz") {
      Order[5] = i;
      ++BDIM;
    } else {
      std::cerr << "ERROR: Incorrect format" << std::endl;
      throw;
    }

    ++i;
  }


  if (XDIM > 3 || BDIM > 3) {
    std::cerr << "ERROR: spatial or B-field dimensions are too large(>3)" << std::endl;
    throw;
  }

  std::vector<size_t> PrefOrder;
  for (size_t j = 0; j != 6; ++j) {
    if (Order[j] != -1) {
      HasXB[i] = true;

      PrefOrder.push_back(Order[j]);
    }
  }


  std::vector< std::vector<double> > D;

  std::vector<double> v(PrefOrder.size(), 0);

  while (!fi.eof()) {
    for (size_t j = 0; j != PrefOrder.size(); ++j) {
      fi >> v[ PrefOrder[j] ];
    }

    if (fi.eof()) {
      break;
    }

    for (size_t k = 0; k != PrefOrder.size(); ++k) {
      std::cout << v[k] << " ";
    }
    std::cout << std::endl;

    D.push_back(v);
  }


  if (XDIM == 1) {
    std::sort(D.begin(), D.end(), SortVector1);
  } else if (XDIM == 2) {
    std::sort(D.begin(), D.end(), SortVector2);
  } else if (XDIM == 3) {
    std::cout << 3 << std::endl;
    std::sort(D.begin(), D.end(), SortVector3);
  } else {
    std::cerr << "ERROR: XDIM incorrect" << std::endl;
    throw;
  }


  for (size_t i = 0; i != D.size(); ++i) {
    for (size_t j = 0; j < XDIM + BDIM; ++j) {
      printf("%15.6E ", D[i][j]);
    }
    std::cout << std::endl;
  }



  // Number of points in each physical dimension for regularized grid
  size_t const NPX =    20;
  size_t const NPY =    20;
  size_t const NPZ =    20;
  size_t const NPTotal = NPX * NPY * NPZ;

  // When you're dealing with potentially a huge amount of memory it's better to do this yourself...
  std::vector<double>* RegularizedData = new std::vector<double>(NPTotal * (XDIM + BDIM), 0);

  std::cout << RegularizedData->size() << std::endl;

  sleep(200);

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [InFile] [format: X Bx Y By]" << std::endl;
    return 1;
  }

  ReadFile(argv[1], argv[2]);

  return 0;
}
