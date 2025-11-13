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

#include "pylith/feassemble/feassemblefwd.hh"

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/problems/problemsfwd.hh" // HASA Physics
#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/problems/Observer.hh" // USES Observer

#include <memory> // ISA std::enable_shared_from_this
#include "pylith/utils/array.hh" // HASA pylith::integer_array
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

class pylith::feassemble::PhysicsImplementation : public pylith::utils::GenericComponent, public std::enable_shared_from_this<pylith::feassemble::PhysicsImplementation> {
    friend class TestPhysicsImplementation; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor.
     *
     * @param[in] physics Physics implemented by constraint.
     */
    PhysicsImplementation(const std::shared_ptr<pylith::problems::Physics>& physics);

    /** Constructor for std::shared_ptr.
     *
     * @param[in] physics Physics implemented by constraint.
     */
    static
    std::shared_ptr<PhysicsImplementation> create(const std::shared_ptr<pylith::problems::Physics>& physics);

    /// Destructor.
    virtual ~PhysicsImplementation(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Get mesh associated with domain governed by specified physics.
     *
     * @returns Mesh associated with constrained domain.
     */
    virtual
    const pylith::topology::Mesh& getPhysicsDomainMesh(void) const = 0;

    /** Get name of label marking material.
     *
     * @returns Name of label for material (from mesh generator).
     */
    const char* getPhysicsLabelName(void) const;

    /** Get value of label marking material.
     *
     * @returns Value of label for material (from mesh generator).
     */
    int getPhysicsLabelValue(void) const;

    /** Get auxiliary field.
     *
     * @returns field Field with auxiliary subfields.
     */
    const pylith::topology::Field* getAuxiliaryField(void) const;

    /** Get diagnostic field.
     *
     * @returns field Field with diagnostic subfields.
     */
    const pylith::topology::Field* getDiagnosticField(void) const;

    /** Get derived field.
     *
     * @return field Field with subfields derived from solution.
     */
    const pylith::topology::Field* getDerivedField(void) const;

    /** Notify observers of current solution.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     */
    void notifyObservers(const pylith::real t,
                         const pylith::integer tindex,
                         const pylith::topology::Field& solution,
                         const pylith::problems::Observer::NotificationType notification);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    std::shared_ptr<pylith::problems::Physics> _physics;
    std::shared_ptr<pylith::problems::ObserversPhysics> _observers;

    std::shared_ptr<pylith::topology::Field> _auxiliaryField;
    std::shared_ptr<pylith::topology::Field> _diagnosticField;
    std::shared_ptr<pylith::topology::Field> _derivedField;

    std::unique_ptr<pylith::utils::EventLogger> _logger;

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    PhysicsImplementation(void) = delete;
    PhysicsImplementation(const PhysicsImplementation &) = delete;
    const PhysicsImplementation& operator=(const PhysicsImplementation&) = delete;

}; // class PhysicsImplementation

// End of file
