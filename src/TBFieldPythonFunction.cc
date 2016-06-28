#include "TBFieldPythonFunction.h"


TBFieldPythonFunction::TBFieldPythonFunction (PyObject* Function)
{
  Py_INCREF(Function);
  fPythonFunction = Function;
}



TBFieldPythonFunction::~TBFieldPythonFunction ()
{
  Py_DECREF(fPythonFunction);
}




TVector3D TBFieldPythonFunction::GetB (TVector3D const& X) const
{
  // Create a python list
  PyObject *InputList = PyList_New(0);
  PyList_Append(InputList, Py_BuildValue("d", X.GetX()));
  PyList_Append(InputList, Py_BuildValue("d", X.GetY()));
  PyList_Append(InputList, Py_BuildValue("d", X.GetZ()));


  if (!PyCallable_Check(fPythonFunction)) {
    std::cout << "PyCallable_Check fail" << std::endl;
    throw;
  }

  PyObject* InputTuple;
  InputTuple = Py_BuildValue("(O)", InputList);

  PyObject* OutputTuple = PyEval_CallObject(fPythonFunction, InputTuple);
  Py_DECREF(InputTuple);

  if (OutputTuple == NULL) {
    std::cout << "No object" << std::endl;
    throw;
  }

  PyObject* OutputList;
  if (!PyArg_Parse(OutputTuple, "O!", &PyList_Type, &OutputList)) {
    std::cout << "Didnt make it" << std::endl;
    throw;
  }



  // Observation point
  TVector3D ReturnVector(PyFloat_AsDouble(PyList_GetItem(OutputList, 0)),
                   PyFloat_AsDouble(PyList_GetItem(OutputList, 1)),
                   PyFloat_AsDouble(PyList_GetItem(OutputList, 2)));

  Py_DECREF(OutputTuple);
  Py_DECREF(OutputList);
  Py_DECREF(InputList);

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
