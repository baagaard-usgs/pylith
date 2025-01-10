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

#include "PhysicsStub.hh" // Implementation of class methods

#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "tests/src/StubMethodTracker.hh" // USES StubMethodTracker
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::PhysicsStub::PhysicsStub(void) :
    _auxiliaryFactory(new pylith::feassemble::AuxiliaryFactory) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::PhysicsStub::~PhysicsStub(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::PhysicsStub::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
    Physics::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::problems::PhysicsStub::verifyConfiguration(const pylith::topology::Field& solution) const {
    pylith::testing::StubMethodTracker tracker("pylith::problems::PhysicsStub::verifyConfiguration");
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::problems::PhysicsStub::createIntegrator(const pylith::topology::Field& solution) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::PhysicsStub::createIntegrator");

    return NULL;
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
std::vector<pylith::feassemble::Constraint*>
pylith::problems::PhysicsStub::createConstraints(const pylith::topology::Field& solution) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::PhysicsStub::createConstraints");
    std::vector<pylith::feassemble::Constraint*> constraintArray;

    PYLITH_METHOD_RETURN(constraintArray);
} // createConstraints


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::problems::PhysicsStub::createAuxiliaryField(const pylith::topology::Field& solution,
                                                    const pylith::topology::Mesh& physicsMesh) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::PhysicsStub::createAuxiliaryField");

    return NULL;
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::problems::PhysicsStub::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


// End of file
