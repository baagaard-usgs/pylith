// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================
// SWIG interface
%module meshio

// Header files for module C++ code
%{
#include "pylith/meshio/MeshIO.hh"
#include "pylith/meshio/MeshIOAscii.hh"
#include "pylith/meshio/MeshIOPetsc.hh"
#if defined(ENABLE_CUBIT)
#include "pylith/meshio/MeshIOCubit.hh"
#endif
#include "pylith/meshio/MeshConverter.hh"

#include "pylith/meshio/OutputTrigger.hh"
#include "pylith/meshio/OutputTriggerStep.hh"
#include "pylith/meshio/OutputTriggerTime.hh"
#include "pylith/meshio/DataWriter.hh"
#include "pylith/meshio/DataWriterVTK.hh"
#if defined(ENABLE_HDF5)
#include "pylith/meshio/DataWriterHDF5.hh"
#include "pylith/meshio/DataWriterHDF5Ext.hh"
#endif
#include "pylith/meshio/OutputObserver.hh"
#include "pylith/meshio/OutputSoln.hh"
#include "pylith/meshio/OutputSolnDomain.hh"
#include "pylith/meshio/OutputSolnBoundary.hh"
#include "pylith/meshio/OutputSolnPoints.hh"
#include "pylith/meshio/OutputPhysics.hh"

#include "pylith/utils/arrayfwd.hh"
%}

%include "exception.i"
%exception {
  try {
    $action
  } catch (const std::exception& err) {
    SWIG_exception(SWIG_RuntimeError, err.what());
  } // try/catch
 } // exception

%include "std_string.i"
%include "typemaps.i"
%include "../include/scalartypemaps.i"
%include "../include/chararray.i"

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
%}
%include "../include/numpy.i"
%init %{
import_array();
%}

// Interfaces
%include "../utils/PyreComponent.i"
%include "../problems/ObserverSoln.i"
%include "../problems/ObserverPhysics.i"
%include "MeshIOObj.i"
%include "MeshIOAscii.i"
%include "MeshIOPetsc.i"
#if defined(ENABLE_CUBIT)
%include "MeshIOCubit.i"
#endif
%include "MeshConverter.i"

%include "OutputTrigger.i"
%include "OutputTriggerStep.i"
%include "OutputTriggerTime.i"
%include "DataWriter.i"
%include "DataWriterVTK.i"
#if defined(ENABLE_HDF5)
%include "DataWriterHDF5.i"
%include "DataWriterHDF5Ext.i"
#endif
%include "OutputObserver.i"
%include "OutputSoln.i"
%include "OutputSolnDomain.i"
%include "OutputSolnBoundary.i"
%include "OutputSolnPoints.i"
%include "OutputPhysics.i"

// End of file
