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

/*
 * We cast the problem in terms of F(t,s,\dot{s}) = G(t,s), s(t0) = s0.
 *
 * In PETSc time stepping (TS) notation, G is the RHS, and F is the I
 * function (which we call the LHS).
 */

#include "pylith/problems/problemsfwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/feassemble/feassemblefwd.hh" // HOLDSA Integrator, Constraint, Observers
#include "pylith/materials/materialsfwd.hh" // HOLDSA Material
#include "pylith/bc/bcfwd.hh" // HOLDSA BoundaryCondition
#include "pylith/faults/faultsfwd.hh" // HOLDSA FaultCohesive
#include "pylith/testing/testingfwd.hh" // MMSTest ISA friend
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HASA GravityField

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field
#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional

#include "pylith/utils/petscfwd.h" // USES PetscVec, PetscMat

#include "pylith/problems/Physics.hh" // USES Problem::Formulation

#include "pylith/utils/array.hh" // HASA std::vector

class pylith::problems::Problem : public pylith::utils::PyreComponent {
    friend class TestProblem; // unit testing
    friend class pylith::testing::MMSTest; // MMS testing

    // PUBLIC ENUM /////////////////////////////////////////////////////////////////////////////////////////////////////
public:

    enum SolverTypeEnum {
        LINEAR, // Linear solver.
        NONLINEAR, // Nonlinear solver.
    }; // SolverType

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    Problem(void);

    /// Destructor
    virtual ~Problem(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set formulation for equations.
     *
     * @param[in] value Formulation type.
     */
    void setFormulation(const pylith::problems::Physics::FormulationEnum value);

    /** Get formulation for equations.
     *
     * @returns Formulation type.
     */
    pylith::problems::Physics::FormulationEnum getFormulation(void) const;

    /** Set solver type.
     *
     * @param[in] value Solver type.
     */
    void setSolverType(const SolverTypeEnum value);

    /** Get solver type.
     *
     * @returns Solver type.
     */
    SolverTypeEnum getSolverType(void) const;

    /** Specify which default PETSc options to use.
     *
     * @param[in] flags Flags indicating which default PETSc options to set.
     */
    void setPetscDefaults(const int flags);

    /** Set manager of scales used to nondimensionalize problem.
     *
     * @param[in] normalizer Nondimensionalizer.
     */
    void setNormalizer(const std::shared_ptr<spatialdata::units::Nondimensional>& normalizer);

    /** Set gravity field.
     *
     * @param[in] gravityField Gravity field.
     */
    void setGravityField(const std::shared_ptr<spatialdata::spatialdb::GravityField>& const g);

    /** Register observer to receive notifications.
     *
     * Observers are used for output.
     *
     * @param[in] observer Observer to receive notifications.
     */
    void registerObserver(const std::shared_ptr<pylith::problems::ObserverSoln>& observer);

    /** Remove observer from receiving notifications.
     *
     * @param[in] observer Observer to remove.
     */
    void removeObserver(const std::shared_ptr<pylith::problems::ObserverSoln>& observer);

    /** Set solution field.
     *
     * @param[in] field Solution field.
     */
    void setSolution(const std::shared_ptr<pylith::topology::Field>& field);

    /** Get solution field.
     *
     * @returns Solution field.
     */
    const pylith::topology::Field* getSolution(void) const;

    /** Get time derivative solution field.
     *
     * @returns Time derivative of solution field.
     */
    const pylith::topology::Field* getSolutionDot(void) const;

    /** Set materials.
     *
     * @param[in] materials Array of materials.
     */
    void setMaterials(const std::vector<std::shared_ptr<pylith::materials::Material> >& materials);

    /** Set boundary conditions.
     *
     * @param[in] bcs Array of boundary conditions.
     */
    void setBoundaryConditions(const std::vector<std::shared_ptr<pylith::bc::BoundaryCondition> >& bcs);

    /** Set interior interface conditions.
     *
     * @param[in] interfaces Array of interior interfaces.
     */
    void setInterfaces(const std::vector<std::shared_ptr<pylith::faults::FaultCohesive> >& interfaces);

    /** Do minimal initialization.
     *
     * @param mesh Finite-element mesh.
     */
    virtual
    void preinitialize(const pylith::topology::Mesh& mesh);

    /// Verify configuration.
    virtual
    void verifyConfiguration(void) const;

    /// Initialize problem.
    virtual
    void initialize(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    std::unique_ptr<pylith::feassemble::IntegrationData> _integrationData; /// > Data needed to integrate PDE.

    std::shared_ptr<spatialdata::units::Nondimensional> _normalizer; ///< Nondimensionalization of scales.
    std::shared_ptr<spatialdata::spatialdb::GravityField> _gravityField; ///< Gravity field.

    std::vector<std::shared_ptr<pylith::materials::Material> > _materials; ///< Array of materials.
    std::vector<std::shared_ptr<pylith::bc::BoundaryCondition> > _bcs; ///< Array of boundary conditions.
    std::vector<std::shared_ptr<pylith::faults::FaultCohesive> > _interfaces; ///< Array of interior interfaces.

    std::vector<std::unique_ptr<pylith::feassemble::Integrator> > _integrators; ///< Array of integrators.
    std::vector<std::unique_ptr<pylith::feassemble::Constraint> > _constraints; ///< Array of constraints.
    std::unique_ptr<pylith::problems::ObserversSoln> _observers; ///< Subscribers of solution updates.

    pylith::problems::Physics::FormulationEnum _formulation; ///< Formulation for equations.
    SolverTypeEnum _solverType; ///< Problem (solver) type.
    int _petscDefaults; ///< Flags for PETSc default options for problem.

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /// Check material and interface label values.
    void _checkMaterialLabels(void) const;

    /// Create array of integrators from materials, interfaces, and boundary conditions.
    void _createIntegrators(void);

    /// Create array of constraints from materials, interfaces, and boundary conditions.
    void _createConstraints(void);

    /// Setup solution subfields and discretization.
    void _setupSolution(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    Problem(const Problem&) = delete;
    const Problem& operator=(const Problem&) = delete;

}; // Problem

// End of file
