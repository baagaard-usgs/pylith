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

#include "ObserverSolnStub.hh" // Implementation of class methods

#include "tests/src/StubMethodTracker.hh" // USES StubMethodTracker

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::ObserverSolnStub::ObserverSolnStub(void) :
    _timeScale(1.0) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ObserverSolnStub::~ObserverSolnStub(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Set time scale.
void
pylith::problems::ObserverSolnStub::setTimeScale(const pylith::real value) {
    _timeScale = value;
} // setTimeScale


// ---------------------------------------------------------------------------------------------------------------------
// Get time scale.
pylith::real
pylith::problems::ObserverSolnStub::getTimeScale(void) const {
    return _timeScale;
} // getTimeScale


// ---------------------------------------------------------------------------------------------------------------------
// Verify observer is compatible with solution.
void
pylith::problems::ObserverSolnStub::verifyConfiguration(const pylith::topology::Field& solution) const {
    pylith::testing::StubMethodTracker tracker("pylith::problems::ObserverPhysicsStub::verifyConfiguration");
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Receive update (subject of observer).
void
pylith::problems::ObserverSolnStub::update(const pylith::real t,
                                           const pylith::integer tindex,
                                           const pylith::topology::Field& solution,
                                           const NotificationType notification) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::ObserverPhysicsStub::update");
} // update


// End of file
