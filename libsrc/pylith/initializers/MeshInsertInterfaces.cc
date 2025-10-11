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

#include "pylith/initializers/MeshInsertInterfaces.hh" // implementation of class methods

#include "pylith/problems/Problem.hh" // USES Problem
#include "pylith/materials/Material.hh" // USES Material
#include "pylith/faults/FaultCohesive.hh" // USES FaultCohesive
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ------------------------------------------------------------------------------------------------
// Default constructor
pylith::initializers::MeshInsertInterfaces::MeshInsertInterfaces(void) {}


// ------------------------------------------------------------------------------------------------
// Default destructor
pylith::initializers::MeshInsertInterfaces::~MeshInsertInterfaces(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::initializers::MeshInsertInterfaces::deallocate(void) {
    InitializePhase::deallocate();
} // deallocate


// ------------------------------------------------------------------------------------------------
// Run initialization phase.
pylith::topology::Mesh*
pylith::initializers::MeshInsertInterfaces::run(pylith::topology::Mesh* mesh,
                                                const pylith::problems::Problem& problem) {
    PYLITH_METHOD_BEGIN;

    PylithInt cohesiveLabelValue = 100;
    for (auto material : problem.getMaterials()) {
        const PylithInt materialLabelValue = material->getLabelValue();
        cohesiveLabelValue = std::max(cohesiveLabelValue, materialLabelValue);
    } // for

    for (auto interface : problem.getInterfaces()) {
        interface->setCohesiveLabelValue(cohesiveLabelValue);
        interface->transformTopology(mesh);
        cohesiveLabelValue += 1;
    } // for

    PYLITH_METHOD_RETURN(mesh);
} // run


// End of file
