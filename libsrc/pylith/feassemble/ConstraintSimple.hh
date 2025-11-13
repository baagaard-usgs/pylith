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

#include "pylith/feassemble/Constraint.hh" // ISA Constraint

#include "pylith/problems/problemsfwd.hh" // HASA Physics
#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/utils/array.hh" // HASA pylith::integer_array
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

class pylith::feassemble::ConstraintSimple : public pylith::feassemble::Constraint {
    // friend class TestConstraintSimple; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Factory for std::shared_ptr.
     *
     * @param[in] physics Physics implemented by constraint.
     */
    static
    std::shared_ptr<ConstraintSimple> create(const std::shared_ptr<pylith::problems::Physics>& physics);

    /// Destructor.
    virtual ~ConstraintSimple(void) override;

    /** Set user function specifying constrained values.
     *
     * @param[in] fn Function specifying contrained values.
     */
    void setUserFn(const PetscUserFieldFunc fn);

    /** Initialize constraint.
     *
     * @param[in] solution Solution field (layout).
     */
    virtual
    void initialize(const pylith::topology::Field& solution) override;

    /** Set constrained values in solution field.
     *
     * @param[inout] integrationData Data needed to integrate governing equation.
     */
    virtual
    void setSolution(pylith::feassemble::IntegrationData& integrationData) override;

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
protected:

    /** Constructor only used by factory.
     *
     * @param[in] physics Physics implemented by constraint.
     */
    ConstraintSimple(const std::shared_ptr<pylith::problems::Physics>& physics);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    PetscUserFieldFunc _fn; ///< Function for computing constrained values.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    ConstraintSimple(void) = delete;
    ConstraintSimple(const Constraint &) = delete;
    const ConstraintSimple& operator=(const ConstraintSimple&) = delete;

}; // class Constraint.

// End of file
