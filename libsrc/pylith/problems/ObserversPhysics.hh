// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/problems/problemsfwd.hh" // forward declarations

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/problems/ObserverPhysics.hh" // USES ObserverPhysics
#include "pylith/feassemble/feassemblefwd.hh" // USES PhysicsImplementation
#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/types.hh" // USES pylith::real, pylith::integer

#include <set> // USES std::set

class pylith::problems::ObserversPhysics : public pylith::utils::GenericComponent {
    friend class TestObserversPhysics; // unit testing
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    ObserversPhysics(void);

    /// Destructor
    virtual ~ObserversPhysics(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Register observer to receive notifications.
     *
     * ObserversPhysics are used for output.
     *
     * @param[in] observer Observer to receive notifications.
     */
    void registerObserver(const std::shared_ptr<pylith::problems::ObserverPhysics>& observer);

    /** Remove observer from receiving notifications.
     *
     * @param[in] observer Observer to remove.
     */
    void removeObserver(const std::shared_ptr<pylith::problems::ObserverPhysics>& observer);

    /** Get number of observers.
     *
     * @returns Number of observers.
     */
    size_t size(void) const;

    /** Set physics implementation in observers (for callbacks)
     *
     * @param[in] physics Physics implementation being observed.
     */
    void setPhysicsImplementation(const std::shared_ptr<pylith::feassemble::PhysicsImplementation>& physics);

    /** Set time scale in observers.
     *
     * @param[in] value Time scale for dimensionalizing time.
     */
    void setTimeScale(const pylith::real value);

    /** Verify observers are compatible.
     *
     * @param[in] solution Solution field.
     */
    void verifyObservers(const pylith::topology::Field& solution) const;

    /** Send observers an update.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     * @param[in] notification Type of notification.
     */
    void notifyObservers(const pylith::real t,
                         const pylith::integer tindex,
                         const pylith::topology::Field& solution,
                         const pylith::problems::Observer::NotificationType notification);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::set<std::shared_ptr<pylith::problems::ObserverPhysics> > _observers; ///< Subscribers of updates.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ObserversPhysics(const ObserversPhysics&) = delete;
    const ObserversPhysics& operator=(const ObserversPhysics&) = delete;

};

// ObserversPhysics

// End of file
