#include "TBFieldPythonFunction.h"

// UPDATE: exceptions

TBFieldPythonFunction::TBFieldPythonFunction (PyObject* Function)
{
  // Constructor takes a python object, which should be a function
  // Increment reference because we're going to keep it..

  Py_INCREF(Function);
  fPythonFunction = Function;

  // Check to see the function is callable
  if (!PyCallable_Check(fPythonFunction)) {
    std::cout << "PyCallable_Check fail" << std::endl;
    throw std::invalid_argument("python function not callable");
  }
}



TBFieldPythonFunction::~TBFieldPythonFunction ()
{
  // When exit, decrement reference since we're done with it
  Py_DECREF(fPythonFunction);
}




TVector3D TBFieldPythonFunction::GetB (TVector3D const& X) const
{
  // Create a python list from input vector
  PyObject *InputList = PyList_New(0);
  PyList_Append(InputList, Py_BuildValue("d", X.GetX()));
  PyList_Append(InputList, Py_BuildValue("d", X.GetY()));
  PyList_Append(InputList, Py_BuildValue("d", X.GetZ()));


  // Check to see the function is callable
  if (!PyCallable_Check(fPythonFunction)) {
    std::cout << "PyCallable_Check fail" << std::endl;
    throw;
  }

  // Build the input object for the python function
  PyObject* InputTuple;
  InputTuple = Py_BuildValue("(O)", InputList);

  // Call python function
  PyObject* OutputTuple = PyEval_CallObject(fPythonFunction, InputTuple);

  // We're done with the input object
  Py_DECREF(InputTuple);

  // If the output is null we didn't get anything
  if (OutputTuple == NULL) {
    std::cout << "No object" << std::endl;
    throw;
  }

  // Get a python list from output tuple
  PyObject* OutputList;
  if (!PyArg_Parse(OutputTuple, "O!", &PyList_Type, &OutputList)) {
    std::cout << "Didnt make it" << std::endl;
    throw;
  }



  // Observation point from python list
  TVector3D ReturnVector(PyFloat_AsDouble(PyList_GetItem(OutputList, 0)),
                         PyFloat_AsDouble(PyList_GetItem(OutputList, 1)),
                         PyFloat_AsDouble(PyList_GetItem(OutputList, 2)));

  // Decrement object references no longer needed
  Py_DECREF(OutputTuple);
  Py_DECREF(OutputList);
  Py_DECREF(InputList);

  // Return the magnetic field vector
  return ReturnVector;
}




TVector3D TBFieldPythonFunction::GetB (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z));
}




double TBFieldPythonFunction::GetBx (double const X, double const Y, double const Z) const
{
  return this->GetB(X, Y, Z).GetX();
}




double TBFieldPythonFunction::GetBy (double const X, double const Y, double const Z) const
{
  return this->GetB(X, Y, Z).GetY();
}




double TBFieldPythonFunction::GetBz (double const X, double const Y, double const Z) const
{
  return this->GetB(X, Y, Z).GetZ();
}
