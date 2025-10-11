// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/initializers/Initializer.hh" // implementation of class methods

#include "pylith/initializers/InitializePhase.hh" // HASA InitializePhase
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*


// ------------------------------------------------------------------------------------------------
// Default constructor
pylith::initializers::Initializer::Initializer(void) {}


// ------------------------------------------------------------------------------------------------
// Default destructor
pylith::initializers::Initializer::~Initializer(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::initializers::Initializer::deallocate(void) {
    _phases.resize(0); // :TODO: Use shared pointer
}


// ------------------------------------------------------------------------------------------------
// Set phases.
void
pylith::initializers::Initializer::setPhases(pylith::initializers::InitializePhase** phases,
                                             const size_t numPhases) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Initializer::setPhases("<<phases<<", numPhases="<<numPhases<<")");

    assert( (!phases && 0 == numPhases) || (phases && 0 < numPhases) );

    _phases.resize(numPhases);
    for (int i = 0; i < numPhases; ++i) {
        _phases[i] = phases[i];
    } // for

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Run initialization phase.
pylith::topology::Mesh*
pylith::initializers::Initializer::initialize(pylith::topology::Mesh* mesh,
                                              const pylith::problems::Problem& problem) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Initializer::initialize(mesh="<<typeid(mesh).name()<<")");
    assert(mesh);

    pylith::topology::Mesh* newMesh = nullptr;
    pylith::topology::Mesh* phaseMesh = mesh;
    const size_t numPhases = _phases.size();
    for (size_t i = 0; i < numPhases; ++i) {
        assert(_phases[i]);

        newMesh = _phases[i]->run(phaseMesh, problem);
        phaseMesh = newMesh;
    } // for

    PYLITH_METHOD_RETURN(newMesh);
} // initialize


// End of file
