////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Jun 17 10:39:12 EDT 2016
//
// This is the python-C extension which allows access to the c++
// class SRSP.
//
////////////////////////////////////////////////////////////////////

// Include Python.h first!
#include <Python.h>

#include "SRS.h"
#include <iostream>
#include <vector>


typedef struct {
  // Define the SRSObject struct which contains the class I want
  PyObject_HEAD
  SRS* obj;
} SRSObject;




static void SRSP_dealloc(SRSObject* self)
{
  // Python needs to know how to deallocate things in the struct

  delete self->obj;
  self->ob_type->tp_free((PyObject*)self);
}




static PyObject* SRSP_new (PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  // Python needs to know how to create things in this struct

  SRSObject* self = (SRSObject*) type->tp_alloc(type, 0);
  if (self != NULL) {

    // Create the new object for self
    self->obj = new SRS();
  }

  // Return myself
  return (PyObject*) self;
}




static PyObject* SRSP_SetCTStart (SRSObject* self, PyObject* arg)
{
  // Set the start time in [m] for caculation

  // Grab the value from input
  double val = PyFloat_AsDouble(arg);

  // Set the object variable
  self->obj->SetCTStart(val);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject* SRSP_GetCTStart (SRSObject* self)
{
  // Get the start time in [m] for calculation

  return Py_BuildValue("d", self->obj->GetCTStart());
}




static PyObject* SRSP_SetCTStop (SRSObject* self, PyObject* arg)
{
  // Set the stop time in [m] for calculations

  // Grab the value
  double val = PyFloat_AsDouble(arg);

  // Set the object variable
  self->obj->SetCTStop(val);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject* SRSP_GetCTStop (SRSObject* self)
{
  // Get the CTStop variable from SRS

  return Py_BuildValue("d", self->obj->GetCTStop());
}




static PyObject* SRSP_SetCTStartStop (SRSObject* self, PyObject* args)
{
  // Set the start and stop times for SRS in [m]

  // Grab the values
  double Start, Stop;
  if (! PyArg_ParseTuple(args, "dd", &Start, &Stop)) {
    return NULL;
  }

  // Set the object variable
  self->obj->SetCTStartStop(Start, Stop);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}






static PyObject* SRSP_GetNPointsTrajectory (SRSObject* self)
{
  // Get the numper of points for trajectory calculaton
  return PyInt_FromSize_t(self->obj->GetNPointsTrajectory());
}




static PyObject* SRSP_SetNPointsTrajectory (SRSObject* self, PyObject* arg)
{
  // Set the number of points for trajectory calculation

  // Grab the value from input
  size_t N = PyInt_AsSsize_t(arg);

  // Set the object variable
  self->obj->SetNPointsTrajectory(N);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}








static PyObject* SRSP_AddMagneticField (SRSObject* self, PyObject* args)
{
  // Set the start and stop times for SRS in [m]

  // Grab the values
  char* FileName;
  char* FileFormat;
  if (! PyArg_ParseTuple(args, "ss", &FileName, &FileFormat)) {
    return NULL;
  }

  // Set the object variable
  self->obj->AddMagneticField(FileName, FileFormat);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}








static PyObject* SRSP_GetBField (SRSObject* self, PyObject* args)
{
  // Set the start and stop times for SRS in [m]

  // Python list object
  PyObject * List;

  // Grab the python list
  if (! PyArg_ParseTuple( args, "O!", &PyList_Type, &List)) {
    return NULL;
  }

  // Has to have the correct number of arguments
  if (PyList_Size(List) != 3) {
    return NULL;
  }


  // Grab the values
  double const X = PyFloat_AsDouble(PyList_GetItem(List, 0));
  double const Y = PyFloat_AsDouble(PyList_GetItem(List, 1));
  double const Z = PyFloat_AsDouble(PyList_GetItem(List, 2));

  // Set the object variable
  TVector3D const B = self->obj->GetB(X, Y, Z);

  // Create a python list
  PyObject *PList = PyList_New(0);

  PyList_Append(PList, Py_BuildValue("f", B.GetX()));
  PyList_Append(PList, Py_BuildValue("f", B.GetY()));
  PyList_Append(PList, Py_BuildValue("f", B.GetZ()));

  // Return the python list
  return PList;
}












static PyGetSetDef SRSP_getseters[] = {
  {"ctstart", (getter)SRSP_GetCTStart, (setter)SRSP_SetCTStart, "Start time for calculation in [m]", NULL},
  {"ctstop",  (getter)SRSP_GetCTStop,  (setter)SRSP_SetCTStop,  "Stop time for calculation in [m]", NULL},
  {"npoints_trajectory",  (getter)SRSP_GetNPointsTrajectory,  (setter)SRSP_SetNPointsTrajectory,  "Total number of points for trajectory", NULL},
  {NULL}  /* Sentinel */
};




static PyMethodDef SRSP_methods[] = {
  // We must tell python about the function we allow access as well as give them nice
  // python names, and tell python the method of input parameters.
  {"set_ctstart",  (PyCFunction)SRSP_SetCTStart,  METH_O,       "set the start time in [m]"},
  {"get_ctstart",  (PyCFunction)SRSP_GetCTStart,  METH_NOARGS,  "get the start time in [m]"},
  {"set_ctstop",   (PyCFunction)SRSP_SetCTStop,   METH_O,       "set the stop time in [m]"},
  {"get_ctstop",   (PyCFunction)SRSP_GetCTStop,   METH_NOARGS,  "get the stop time in [m]"},
  {"get_ctstop",   (PyCFunction)SRSP_GetCTStop,   METH_NOARGS,  "get the stop time in [m]"},
  {"set_npoints_trajectory",   (PyCFunction)SRSP_SetNPointsTrajectory,   METH_O,       "set the total number of points for the trajectory"},
  {"get_npoints_trajectory",   (PyCFunction)SRSP_GetNPointsTrajectory,   METH_NOARGS,  "get the total number of points for the trajectory"},
  {"set_ctstartstop",   (PyCFunction)SRSP_SetCTStartStop,  METH_VARARGS,       "set the start and stop time in [m]"},

  {"add_magnetic_field",   (PyCFunction)SRSP_AddMagneticField,  METH_VARARGS,       "add a magnetic field from a file"},
  {"get_bfield",  (PyCFunction)SRSP_GetBField,  METH_VARARGS,  "get the magnetic field at a given position in space (and someday time?)"},

  {NULL}  /* Sentinel */
};




static PyTypeObject SRSType = {
  // The python object.  Fully defined elsewhere.  only put here what you need,
  // otherwise default values

  PyObject_HEAD_INIT(NULL)
  0,                                        /* ob_size */
  "mod.SRS",                                 /* tp_name */
  sizeof(SRSObject),                         /* tp_basicsize */
  0,                                        /* tp_itemsize */
  (destructor) SRSP_dealloc,                 /* tp_dealloc */
  0,                                        /* tp_print */
  0,                                        /* tp_getattr */
  0,                                        /* tp_setattr */
  0,                                        /* tp_compare */
  0,                                        /* tp_repr */
  0,                                        /* tp_as_number */
  0,                                        /* tp_as_sequence */
  0,                                        /* tp_as_mapping */
  0,                                        /* tp_hash */
  0,                                        /* tp_call */
  0,                                        /* tp_str */
  0,                                        /* tp_getattro */
  0,                                        /* tp_setattro */
  0,                                        /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
  "SRSP class",                              /* tp_doc */
  0,                                        /* tp_traverse */
  0,                                        /* tp_clear */
  0,                                        /* tp_richcompare */
  0,                                        /* tp_weaklistoffset */
  0,                                        /* tp_iter */
  0,                                        /* tp_iternext */
  SRSP_methods,                             /* tp_methods */
  0,                                        /* tp_members */
  SRSP_getseters,                           /* tp_getset */
  0,                                        /* tp_base */
  0,                                        /* tp_dict */
  0,                                        /* tp_descr_get */
  0,                                        /* tp_descr_set */
  0,                                        /* tp_dictoffset */
  0,                                        /* tp_init */
  0,                                        /* tp_alloc */
  SRSP_new,                                  /* tp_new */
};




static PyMethodDef module_methods[] = {
  // I do not need
  {NULL}  /* Sentinel */
};




PyMODINIT_FUNC initSRSP()
{
  // Initialization of the module.


  PyObject* m;

  if (PyType_Ready(&SRSType) < 0)
    return;

  m = Py_InitModule3("SRSP", module_methods, "Example module that creates an extension type.");
  if (m == NULL)
    return;

  Py_INCREF(&SRSType);
  PyModule_AddObject(m, "SRSP", (PyObject *)&SRSType);
}


