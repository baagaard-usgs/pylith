// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/problems/Physics.hh" // implementation of class methods

#include "pylith/topology/SubfieldFactory.hh" // USES SubfieldFactory
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_JMETHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> // USES std::runtime_error

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::Physics::Physics(void) :
    _observers(new pylith::problems::ObserversPhysics),
    _labelName(pylith::topology::Mesh::cells_label_name),
    _labelValue(1) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::Physics::~Physics(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Physics::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    _subfieldFactory.reset();
    _normalizer.reset();
    _observers.reset();
    _kernelConstants.resize(0);

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set name of label marking material.
void
pylith::problems::Physics::setLabelName(const char* value) {
    PYLITH_COMPONENT_DEBUG("setLabelName(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for material label.");
    } // if

    _labelName = value;
} // setLabelName


// ------------------------------------------------------------------------------------------------
// Get name of label marking material.
const char*
pylith::problems::Physics::getLabelName(void) const {
    return _labelName.c_str();
} // getLabelName


// ------------------------------------------------------------------------------------------------
// Set value of label marking material.
void
pylith::problems::Physics::setLabelValue(const int value) {
    _labelValue = value;
} // setLabelValue


// ------------------------------------------------------------------------------------------------
// Get value of label marking material.
int
pylith::problems::Physics::getLabelValue(void) const {
    return _labelValue;
} // getLabelValue


// ------------------------------------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::problems::Physics::setNormalizer(const std::shared_ptr<spatialdata::units::Nondimensional>& normalizer) {
    PYLITH_COMPONENT_DEBUG("setNormalizer(dim="<<typeid(normalizer).name()<<")");

    _normalizer = normalizer;
} // setNormalizer


// ------------------------------------------------------------------------------------------------
/** Get manager of scales used to nondimensionalize problem.
 *
 * @param dim Nondimensionalizer.
 */
const spatialdata::units::Nondimensional&
pylith::problems::Physics::getNormalizer(void) const {
    assert(_normalizer);
    return *_normalizer;
} // getNormalizer


// ------------------------------------------------------------------------------------------------
// Set formulation for equations.
void
pylith::problems::Physics::setFormulation(const FormulationEnum value) {
    _formulation = value;
} // setFormulation


// ------------------------------------------------------------------------------------------------
// Set spatial database for populating auxiliary field.
void
pylith::problems::Physics::setAuxiliaryFieldDB(const std::shared_ptr<spatialdata::spatialdb::SpatialDB>& db) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setAuxiliaryFieldDB(db="<<db<<")");

    _auxiliaryFieldDB = db;

    PYLITH_METHOD_END;
} // setAuxiliaryFieldDB


// ------------------------------------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::problems::Physics::setSubfieldDiscretization(const char* subfieldName,
                                                     const pylith::topology::FieldBase::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setSubfieldDiscretization(subfieldName="<<subfieldName<<")");

    pylith::topology::SubfieldFactory* factory = _getSubfieldFactory();assert(factory);
    factory->setDiscretization(subfieldName, discretization);

    PYLITH_METHOD_END;
} // setSubfieldDiscretization


// ----------------------------------------------------------------------
// Register observer to receive notifications.
void
pylith::problems::Physics::registerObserver(const std::shared_ptr<pylith::problems::ObserverPhysics>& observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("registerObserver(observer="<<typeid(observer).name()<<")");

    assert(_observers);
    _observers->registerObserver(observer);

    PYLITH_METHOD_END;
} // registerObserver


// ----------------------------------------------------------------------
// Remove observer from receiving notifications.
void
pylith::problems::Physics::removeObserver(const std::shared_ptr<pylith::problems::ObserverPhysics>& observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("removeObserver(observer="<<typeid(observer).name()<<")");

    assert(_observers);
    _observers->removeObserver(observer);

    PYLITH_METHOD_END;
} // removeObserver


// ------------------------------------------------------------------------------------------------
// Get observers receiving notifications of physics updates.
const std::shared_ptr<pylith::problems::ObserversPhysics>&
pylith::problems::Physics::getObservers(void)  const {
    return _observers;
} // getObservers


// ------------------------------------------------------------------------------------------------
// Get constants used in kernels (point-wise functions).
const pylith::real_array&
pylith::problems::Physics::getKernelConstants(const pylith::real dt) {
    _updateKernelConstants(dt);

    return _kernelConstants;
} // getKernelConstants


// ------------------------------------------------------------------------------------------------
// Create diagnostic field.
std::shared_ptr<pylith::topology::Field>
pylith::problems::Physics::createDiagnosticField(const pylith::topology::Field& solution,
                                                 const pylith::topology::Mesh& physicsMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDiagnosticField(solution="<<solution.getName()<<", physicsMesh=)"<<typeid(physicsMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(nullptr);

}


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
std::shared_ptr<pylith::topology::Field>
pylith::problems::Physics::createDerivedField(const pylith::topology::Field& solution,
                                              const pylith::topology::Mesh& physicsMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getName()<<", physicsMesh=)"<<typeid(physicsMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(nullptr);
} // createDerivedField


// ------------------------------------------------------------------------------------------------
// Update auxiliary field for given time.
void
pylith::problems::Physics::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateAuxiliaryField(auxiliaryField="<<auxiliaryField<<", t="<<t<<") empty method");

    PYLITH_METHOD_END;
} // updateAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::problems::Physics::_updateKernelConstants(const pylith::real dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_updateKernelConstants(dt="<<dt<<")");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // _updateKernelConstants


// End of file
