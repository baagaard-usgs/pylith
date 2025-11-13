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

class pylith::feassemble::ConstraintSpatialDB : public pylith::feassemble::Constraint {
    friend class TestConstraintSpatialDB; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Factory for std::shared_ptr.
     *
     * @param[in] physics Physics implemented by constraint.
     */
    static
    std::shared_ptr<ConstraintSpatialDB> create(const std::shared_ptr<pylith::problems::Physics>& physics);

    /// Destructor.
    ~ConstraintSpatialDB(void) override;

    /** Set constraint kernel.
     *
     * @param kernel Kernel to compute constrained value from auxiliary field.
     */
    void setKernelConstraint(const PetscBdPointFunc kernel);

    /** Initialize constraint.
     *
     * @param[in] solution Solution field (layout).
     */
    void initialize(const pylith::topology::Field& solution) override;

    /** Set auxiliary field values for current time.
     *
     * @param[in] t Current time.
     */
    void setState(const pylith::real t) override;

    /** Set constrained values in solution field.
     *
     * @param[inout] integrationData Data needed to integrate governing equation.
     */
    void setSolution(pylith::feassemble::IntegrationData& integrationData) override;

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Constructor only used by factory.
     *
     * @param[in] physics Physics implemented by constraint.
     */
    ConstraintSpatialDB(const std::shared_ptr<pylith::problems::Physics>& physics);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    PetscBdPointFunc _kernelConstraint; ///< Kernel for computing constrained values from auxiliary field.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    ConstraintSpatialDB(void) = delete;
    ConstraintSpatialDB(const ConstraintSpatialDB &) = delete;
    const ConstraintSpatialDB& operator=(const ConstraintSpatialDB&) = delete;

}; // class Constraint

// End of file
