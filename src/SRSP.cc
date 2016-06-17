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
  // Set only one variable in the class object within the struct.
  // Parsing is different from generic case.  METH_O must be specified

  // Grab the int
  double val = PyFloat_AsDouble(arg);

  // Set the object variable
  self->obj->SetCTStart(val);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




static PyObject* SRSP_GetCTStart (SRSObject* self)
{
  // Get a variable from the class object within the struct and return
  // it as a python object

  return Py_BuildValue("f", self->obj->GetCTStart());
}




static PyMethodDef SRSP_methods[] = {
  // We must tell python about the function we allow access as well as give them nice
  // python names, and tell python the method of input parameters.
  {"setA",  (PyCFunction)SRSP_SetCTStart,  METH_O,       "set the start time in [m]"},
  {"getA",  (PyCFunction)SRSP_GetCTStart,  METH_NOARGS,  "get the start time in [m]"},
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
  SRSP_methods,                              /* tp_methods */
  0,                                        /* tp_members */
  0,                                        /* tp_getset */
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


