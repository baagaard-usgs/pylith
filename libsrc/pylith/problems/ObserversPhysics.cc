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

#include "pylith/problems/ObserversPhysics.hh" // Implementation of class methods

#include "pylith/feassemble/PhysicsImplementation.hh" // USES PhysicsImplementation
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_DEBUG

#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::ObserversPhysics::ObserversPhysics(void) {
    // GenericComponent::setName("observersphysics");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ObserversPhysics::~ObserversPhysics(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::ObserversPhysics::deallocate(void) {
    _observers.clear();
} // deallocate


// ------------------------------------------------------------------------------------------------
// Register observer to receive notifications.
void
pylith::problems::ObserversPhysics::registerObserver(const std::shared_ptr<pylith::problems::ObserverPhysics>& observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("registerObserver(observer="<<typeid(observer).name()<<")");

    if (observer) {
        _observers.insert(observer);
    } // if

    PYLITH_METHOD_END;
} // registerObserver


// ------------------------------------------------------------------------------------------------
// Remove observer from receiving notifications.
void
pylith::problems::ObserversPhysics::removeObserver(const std::shared_ptr<pylith::problems::ObserverPhysics>& observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("removeObserver(observer="<<typeid(observer).name()<<")");

    if (observer) {
        _observers.erase(observer);
    } // if

    PYLITH_METHOD_END;
} // removeObserver


// ------------------------------------------------------------------------------------------------
// Get number of observers.
size_t
pylith::problems::ObserversPhysics::size(void) const {
    return _observers.size();
} // count`


// ------------------------------------------------------------------------------------------------
// Set physics implementation in observers (for callbacks)
void
pylith::problems::ObserversPhysics::setPhysicsImplementation(const std::shared_ptr<pylith::feassemble::PhysicsImplementation>& physics) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setPhysicsImplementation(physics="<<physics<<")");

    for (auto observer : _observers) {
        assert(observer);
        observer->setPhysicsImplementation(physics);
    } // for

    PYLITH_METHOD_END;
} // setPhysicsImplemetation


// ------------------------------------------------------------------------------------------------
// Set time scale in observers.
void
pylith::problems::ObserversPhysics::setTimeScale(const pylith::real value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setTimeScale(value="<<value<<")");

    for (auto observer : _observers) {
        assert(observer);
        observer->setTimeScale(value);
    } // for

    PYLITH_METHOD_END;
} // setTimeScale


// ------------------------------------------------------------------------------------------------
// Verify observers.
void
pylith::problems::ObserversPhysics::verifyObservers(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("verifyObservers(solution="<<solution.getName()<<")");

    for (auto observer : _observers) {
        assert(observer);
        observer->verifyConfiguration(solution);
    } // for

    PYLITH_METHOD_END;
} // verifyObservers


// ------------------------------------------------------------------------------------------------
// Notify observers.
void
pylith::problems::ObserversPhysics::notifyObservers(const pylith::real t,
                                                    const pylith::integer tindex,
                                                    const pylith::topology::Field& solution,
                                                    const pylith::problems::Observer::NotificationType notification) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("notifyObservers(t="<<t<<", tindex="<<tindex<<", solution="<<solution.getName()<<", " <<notification<<")");

    for (auto observer : _observers) {
        assert(observer);
        observer->update(t, tindex, solution, notification);
    } // for

    PYLITH_METHOD_END;
} // notifyObservers


// End of file
