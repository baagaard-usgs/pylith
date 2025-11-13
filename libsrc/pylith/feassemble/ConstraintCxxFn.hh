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

class pylith::feassemble::ConstraintCxxFn : public pylith::feassemble::Constraint {
    friend class TestConstraintCxxFn; // unit testing
    friend class _ConstraintCxxFn; /// private utility class

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Factory for std::shared_ptr.
     *
     * @param[in] physics Physics implemented by constraint.
     */
    static
    std::shared_ptr<ConstraintCxxFn> create(const std::shared_ptr<pylith::problems::Physics>& physics);

    /// Destructor.
    virtual ~ConstraintCxxFn(void) override;

    /** Set user function specifying constrained values.
     *
     * @param[in] fn Function specifying contrained values.
     */
    void setCxxFn(const PetscUserFieldFunc fn);

    /** Set user function time derivative specifying constrained values.
     *
     * @param[in] fnDot Function specifying contrained values time derivative.
     */
    void setCxxFnDot(const PetscUserFieldFunc fnDot);

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
    void setSolution(pylith::feassemble::IntegrationData& integrationData) override;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    /** Constructor
     *
     * @param[in] physics Physics implemented by constraint.
     */
    ConstraintCxxFn(const std::shared_ptr<pylith::problems::Physics>& physics);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    PetscUserFieldFunc _fn; ///< Function for computing constrained values.
    PetscUserFieldFunc _fnDot; ///< Function for computing constrained values.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    ConstraintCxxFn(void) = delete;
    ConstraintCxxFn(const Constraint &) = delete;
    const ConstraintCxxFn& operator=(const ConstraintCxxFn&) = delete;

}; // class ConstraintCxxFn

// End of file
