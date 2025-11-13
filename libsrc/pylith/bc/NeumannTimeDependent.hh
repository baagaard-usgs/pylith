// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
//
#pragma once

/*
 * @brief C++ implementation of Neumann (e.g., traction) boundary condition
 * with time-dependent expression for value.
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

class pylith::bc::NeumannTimeDependent : public pylith::bc::BoundaryCondition {
    friend class TestNeumannTimeDependent; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Factory for std::shared_ptr.
     *
     * @param[in] physics Physics implemented by constraint.
     */
    static
    std::shared_ptr<NeumannTimeDependent> create(void);

    /// Destructor.
    ~NeumannTimeDependent(void) override;

    /// Deallocate PETSc and local data structures.
    void deallocate(void) override;

    /** Set time history database.
     *
     * @param[in] db Time history database.
     */
    void setTimeHistoryDB(const std::shared_ptr<spatialdata::spatialdb::TimeHistory>& th);

    /** Get time history database.
     *
     * @preturns Time history database.
     */
    const std::shared_ptr<spatialdata::spatialdb::TimeHistory>& getTimeHistoryDB(void);

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

    /** Set name of scale associated with Neumann boundary
     * condition (e.g., 'pressure' for elasticity).
     *
     * A Neumann boundary condition constrains the gradient in
     * a solution subfield. In some cases the constraint is
     * actually on a scaled version of the gradient as is the
     * case of a Neumann boundary condition for elasticity
     * that constrains boundary tractions.
     *
     * @param value Name of scale for nondimensionalizing Neumann boundary condition.
     */
    void setScaleName(const char* value);

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

    /** Update auxiliary subfields at beginning of time step.
     *
     * @param[out] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const double t) override;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    /// Default constructor.
    NeumannTimeDependent(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    std::shared_ptr<spatialdata::spatialdb::TimeHistory> _dbTimeHistory; ///< Time history database.
    std::string _scaleName; ///< Name of scale associated with Neumann boundary condition.

    bool _useInitial; ///< Use initial value term.
    bool _useRate; ///< Use rate term.
    bool _useTimeHistory; ///< Use time history term.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    NeumannTimeDependent(const NeumannTimeDependent&) = delete;
    const NeumannTimeDependent& operator=(const NeumannTimeDependent&) = delete;

}; // class NeumannTimeDependent

// End of file
