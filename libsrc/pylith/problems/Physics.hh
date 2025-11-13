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

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/Field.hh" // USES FieldBase
#include "pylith/feassemble/feassemblefwd.hh" // USES Integrator, Constraint

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB
#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional

#include <memory> // HASA std::shared_ptr

class pylith::problems::Physics : public pylith::utils::PyreComponent, public std::enable_shared_from_this<Physics> {
    friend class TestPhysics; // unit testing

    // PUBLIC ENUM ////////////////////////////////////////////////////////////////////////////////
public:

    enum FormulationEnum {
        QUASISTATIC, // Without inertia; implicit time stepping.
        DYNAMIC, // With inertia; explicit time stepping).
        DYNAMIC_IMEX, // With inertia; implicit+explicit time stepping).
    }; // FormulationEnum

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    Physics(void);

    /// Destructor
    virtual ~Physics(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set name of label marking material.
     *
     * @param[in] value Name of label for material (from mesh generator).
     */
    void setLabelName(const char* value);

    /** Get name of label marking material.
     *
     * @returns Name of label for material (from mesh generator).
     */
    const char* getLabelName(void) const;

    /** Set value of label marking material.
     *
     * @param[in] value Value of label for material (from mesh generator).
     */
    void setLabelValue(const int value);

    /** Get value of label marking material.
     *
     * @returns Value of label for material (from mesh generator).
     */
    int getLabelValue(void) const;

    /** Set manager of scales used to nondimensionalize problem.
     *
     * @param normalizer Nondimensionalizer.
     */
    void setNormalizer(const std::shared_ptr<spatialdata::units::Nondimensional>& normalizer);

    /** Get manager of scales used to nondimensionalize problem.
     *
     * @returns Nondimensionalizer.
     */
    const spatialdata::units::Nondimensional& getNormalizer(void) const;

    /** Set formulation for equations.
     *
     * @param[in] value Formulation for equations.
     */
    void setFormulation(const FormulationEnum value);

    /** Set spatial database for populating auxiliary field.
     *
     * @param[in] db Spatial database with initial values for auxiliary field.
     */
    void setAuxiliaryFieldDB(const std::shared_ptr<spatialdata::spatialdb::SpatialDB>& db);

    /** Set discretization information for subfield.
     *
     * @param[in] subfieldName Name of subfield.
     * @param[in] discretization Discretization for subfield.
     */
    void setSubfieldDiscretization(const char* subfieldName,
                                   const pylith::topology::FieldBase::Discretization& discretization);

    /** Register observer to receive notifications.
     *
     * Observers are used for output.
     *
     * @param[in] observer Observer to receive notifications.
     */
    void registerObserver(const std::shared_ptr<pylith::problems::ObserverPhysics>& observer);

    /** Remove observer from receiving notifications.
     *
     * @param[in] observer Observer to remove.
     */
    void removeObserver(const std::shared_ptr<pylith::problems::ObserverPhysics>& observer);

    /** Get observers receiving notifications of physics updates.
     *
     * @returns Observers receiving notifications.
     */
    const std::shared_ptr<pylith::problems::ObserversPhysics>& getObservers(void) const;

    /** Get constants used in kernels (point-wise functions).
     *
     * @param[in] solution Solution field.
     * @param[in] dt Current time step.
     *
     * @return Array of constants.
     */
    const pylith::real_array& getKernelConstants(const pylith::real dt);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution) const = 0;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Integrator if applicable, otherwise nullptr.
     */
    virtual
    std::shared_ptr<pylith::feassemble::Integrator> createIntegrator(const pylith::topology::Field& solution) = 0;

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraints if applicable, otherwise nullptr.
     */
    virtual
    std::vector<std::shared_ptr<pylith::feassemble::Constraint> > createConstraints(const pylith::topology::Field& solution) = 0;

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in] physicsMesh Finite-element mesh associated with physics.
     *
     * @returns Auxiliary field if applicable, otherwise nullptr.
     */
    virtual
    std::shared_ptr<pylith::topology::Field> createAuxiliaryField(const pylith::topology::Field& solution,
                                                                  const pylith::topology::Mesh& physicsMesh) = 0;

    /** Create diagnostic field.
     *
     * @param[in] solution Solution field.
     * @param[in] physicsMesh Finite-element mesh associated with physics.
     *
     * @returns Diagnostic field if applicable, otherwise nullptr.
     */
    virtual
    std::shared_ptr<pylith::topology::Field> createDiagnosticField(const pylith::topology::Field& solution,
                                                                   const pylith::topology::Mesh& physicsMesh);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in] physicsMesh Finite-element mesh associated with physics.
     *
     * @returns Derived field if applicable, otherwise nullptr.
     */
    virtual
    std::shared_ptr<pylith::topology::Field> createDerivedField(const pylith::topology::Field& solution,
                                                                const pylith::topology::Mesh& physicsMesh);

    /** Update time-dependent auxiliary field.
     *
     * @param[inout] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    virtual
    void updateAuxiliaryField(pylith::topology::Field* const auxiliaryField,
                              const double t);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    virtual
    pylith::topology::SubfieldFactory* _getSubfieldFactory(void) = 0;

    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    virtual
    void _updateKernelConstants(const pylith::real dt);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    std::unique_ptr<pylith::topology::SubfieldFactory> _subfieldFactory; ///< Factory for creating subfields.
    std::shared_ptr<spatialdata::spatialdb::SpatialDB> _auxiliaryFieldDB; ///< Database for populating auxiliary field.
    std::shared_ptr<spatialdata::units::Nondimensional> _normalizer; ///< Nondimensionalizer.
    FormulationEnum _formulation; ///< Formulation for equations.
    pylith::real_array _kernelConstants; ///< Constants used in finite-element kernels (point-wise functions).

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    std::shared_ptr<pylith::problems::ObserversPhysics> _observers; ///< Subscribers of updates.
    std::string _labelName; ///< Name of label in mesh associated with physics.
    pylith::integer _labelValue; ///< Value of label in mesh associated with physics.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    Physics(const Physics&) = delete;
    const Physics& operator=(const Physics&) = delete;

}; // Physics

// End of file
