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

#include "pylith/bc/BoundaryCondition.hh" // ISA BoundaryCondition

#include "pylith/topology/topologyfwd.hh" // USES Field

class pylith::bc::NeumannCxxFn : public pylith::bc::BoundaryCondition {
    friend class TestNeumannCxxFn; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Factory for std::shared_ptr.
     *
     * @param[in] physics Physics implemented by constraint.
     */
    static
    std::shared_ptr<NeumannCxxFn> create(void);

    /// Destructor.
    ~NeumannCxxFn(void) override;

    /// Deallocate PETSc and local data structures.
    void deallocate(void) override;

    /** Set user function specifying field on boundary.
     *
     * @param[in] fn Function specifying field on boundary.
     */
    void setCxxFn(PetscBdPointFunc fn);

    /** Get user function specifying field on boundary
     *
     * @preturns Function specifying field on boundary.
     */
    PetscBdPointFunc getCxxFn(void) const;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Integrator if applicable, otherwise nullptr.
     */
    std::shared_ptr<pylith::feassemble::Integrator> createIntegrator(const pylith::topology::Field& solution) override;

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise nullptr.
     */
    std::vector<std::shared_ptr<pylith::feassemble::Constraint> > createConstraints(const pylith::topology::Field& solution) override;

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Auxiliary field if applicable, otherwise nullptr.
     */
    std::shared_ptr<pylith::topology::Field> createAuxiliaryField(const pylith::topology::Field& solution,
                                                                  const pylith::topology::Mesh& domainMesh) override;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    /// Default constructor.
    NeumannCxxFn(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    PetscBdPointFunc _fn; ///< Kernel for boundary integration.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    NeumannCxxFn(const NeumannCxxFn&) = delete;
    const NeumannCxxFn& operator=(const NeumannCxxFn&) = delete;

}; // class NeumannCxxFn

// End of file
