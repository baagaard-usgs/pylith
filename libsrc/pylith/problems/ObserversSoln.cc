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

#include "pylith/problems/ObserversSoln.hh" // Implementation of class methods

#include "pylith/problems/ObserverSoln.hh" // USES ObserverSoln
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_DEBUG

#include <typeinfo> // USES typeid()

// ----------------------------------------------------------------------
// Constructor.
pylith::problems::ObserversSoln::ObserversSoln(void) {
    GenericComponent::setName("observerssoln");
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::problems::ObserversSoln::~ObserversSoln(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::ObserversSoln::deallocate(void) {
    _observers.clear();
} // deallocate


// ----------------------------------------------------------------------
// Register observer to receive notifications.
void
pylith::problems::ObserversSoln::registerObserver(const std::shared_ptr<pylith::problems::ObserverSoln>& observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("registerObserver(observer="<<typeid(observer).name()<<")");

    if (observer) {
        observer->index = _observers.size();
        _observers.insert(observer);
    } // if

    PYLITH_METHOD_END;
} // registerObserver


// ----------------------------------------------------------------------
// Remove observer from receiving notifications.
void
pylith::problems::ObserversSoln::removeObserver(const std::shared_ptr<pylith::problems::ObserverSoln>& observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("removeObserver(observer="<<typeid(observer).name()<<")");

    if (observer) {
        _observers.erase(observer);
    } // if

    PYLITH_METHOD_END;
} // removeObserver


// ----------------------------------------------------------------------
// Set time scale in observers.
void
pylith::problems::ObserversSoln::setTimeScale(const pylith::real value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setTimeScale(value="<<value<<")");

    for (auto observer : _observers) {
        assert(observer);
        observer->setTimeScale(value);
    } // for

    PYLITH_METHOD_END;
} // setTimeScale


// ----------------------------------------------------------------------
// Verify observers.
void
pylith::problems::ObserversSoln::verifyObservers(const pylith::topology::Field& solution) const {
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
pylith::problems::ObserversSoln::notifyObservers(const pylith::real t,
                                                 const pylith::integer tindex,
                                                 const pylith::topology::Field& solution,
                                                 const pylith::problems::Observer::NotificationType notification) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("notifyObservers(t="<<t<<", tindex="<<tindex<<", solution="<<solution.getName()<<")");

    for (auto observer : _observers) {
        assert(observer);
        observer->update(t, tindex, solution, notification);
    } // for

    PYLITH_METHOD_END;
} // notifyObservers


// ------------------------------------------------------------------------------------------------
// Comparison function for keeping set of observers in order.
bool
pylith::problems::ObserversSoln::_compare::operator()(const std::shared_ptr<ObserverSoln>& a,
                                                      const std::shared_ptr<ObserverSoln>& b) const {
    assert(a);
    assert(b);
    return a->index < b->index;
} // _compare


// End of file
