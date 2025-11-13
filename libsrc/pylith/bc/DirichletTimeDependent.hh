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
 * @brief C++ implementation of Dirichlet (prescribed values at
 * degrees of freedom) boundary condition with time-dependent expression
 * for value.
 *
 * f(x,t) = f_0(x) + \dot{f}_1(x)*(t-t_1(x)) + f_2(x)*a(t-t_2(x)).
 *
 * Auxiliary fields:
 *     if _useInitial
 *         initial amplitude (scalar or vector) f_0(x)
 *    if _useRate
 *        rate amplitude (scalar or vector) \dot{f}_1(x)
 *        rate start (scalar) t_1(x)
 *    if _useTimeHistory
 *        time history amplitude (scalar or vector) f_2(x)
 *        time history start (scalar) t_2(x)
 *        time history value (scalar) a(t-t_2(x))
 */

#include "pylith/bc/BoundaryCondition.hh" // ISA BoundaryCondition

#include "pylith/topology/topologyfwd.hh" // USES Field

class pylith::bc::DirichletTimeDependent : public pylith::bc::BoundaryCondition {
    friend class TestDirichletTimeDependent; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Factory for std::shared_ptr.
     *
     * @param[in] physics Physics implemented by constraint.
     */
    static
    std::shared_ptr<DirichletTimeDependent> create(void);

    /// Destructor.
    ~DirichletTimeDependent(void) override;

    /// Deallocate PETSc and local data structures.
    void deallocate(void) override;

    /** Set indices of constrained degrees of freedom at each location.
     *
     * Example: [0, 1] to apply forces to x and y degrees of freedom in
     * a Cartesian coordinate system.
     *
     * @param[in] dof Array of indices for constrained degrees of freedom.
     */
    void setConstrainedDOF(const pylith::integer_array& dof);

    /** Get indices of constrained degrees of freedom.
     *
     * @returns Array of indices for constrained degrees of freedom.
     */
    const pylith::integer_array& getConstrainedDOF(void) const;

    /** Set time history database.
     *
     * @param[in] db Time history database.
     */
    void setTimeHistoryDB(const std::shared_ptr<spatialdata::spatialdb::TimeHistory>& th);

    /** Get time history database.
     *
     * @preturns Time history database.
     */
    const spatialdata::spatialdb::TimeHistory& getTimeHistoryDB(void);

    /** Use initial value term in time history expression.
     *
     * @param[in] value True if using initial value term in expression.
     */
    void useInitial(const bool value);

    /** Get flag associated with using initial value term in time history expression.
     *
     * @returns True if using initial value term in expression, false otherwise.
     */
    bool useInitial(void) const;

    /** Use rate value term in time history expression.
     *
     * @param[in] value True if using rate value term in expression.
     */
    void useRate(const bool value);

    /** Get flag associated with using rate value term in time history expression.
     *
     * @returns True if using rate value term in expression, false otherwise.
     */
    bool useRate(void) const;

    /** Use time history term in time history expression.
     *
     * @param[in] value True if using time history term in expression.
     */
    void useTimeHistory(const bool value);

    /** Get flag associated with using time history term in time history expression.
     *
     * @returns True if using time history term in expression, false otherwise.
     */
    bool useTimeHistory(void) const;

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

    /** Update time-dependent auxiliary field.
     *
     * @param[out] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const double t) override;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    /// Default constructor.
    DirichletTimeDependent(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::integer_array _constrainedDOF; ///< List of constrained degrees of freedom at each location.
    std::shared_ptr<spatialdata::spatialdb::TimeHistory> _dbTimeHistory; ///< Time history database.

    bool _useInitial; ///< Use initial value term.
    bool _useRate; ///< Use rate term.
    bool _useTimeHistory; ///< Use time history term.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    DirichletTimeDependent(const DirichletTimeDependent&) = delete;
    const DirichletTimeDependent& operator=(const DirichletTimeDependent&) = delete;

}; // class DirichletTimeDependent

// End of file
