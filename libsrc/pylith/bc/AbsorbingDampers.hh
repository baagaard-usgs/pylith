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

#include "pylith/bc/BoundaryCondition.hh" // ISA Physics

#include "pylith/topology/topologyfwd.hh" // USES Field

class pylith::bc::AbsorbingDampers : public pylith::bc::BoundaryCondition {
    friend class TestAbsorbingDampers; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Factory for std::shared_ptr.
     *
     * @param[in] physics Physics implemented by constraint.
     */
    static
    std::shared_ptr<AbsorbingDampers> create(void);

    /// Destructor.
    ~AbsorbingDampers(void) override;

    /// Deallocate PETSc and local data structures.
    void deallocate(void) override;

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const override;

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

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Default constructor.
    AbsorbingDampers(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    AbsorbingDampers(const AbsorbingDampers&) = delete;
    const AbsorbingDampers& operator=(const AbsorbingDampers&) = delete;

};

// class AbsorbingDampers

// End of file
