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

#include "pylith/bc/DirichletTimeDependent.hh" // implementation of object methods

#include "pylith/bc/SubfieldFactory.hh" // USES SubfieldFactory
#include "pylith/bc/TimeDependentOps.hh" // USES TimeDependentOps

#include "pylith/feassemble/ConstraintSpatialDB.hh" // USES ConstraintSoatialDB
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/fekernels/TimeDependentFn.hh" // USES TimeDependentFn kernels

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace bc {
        class DirichletTimeDependentWrap : public DirichletTimeDependent {
public:

            DirichletTimeDependentWrap(void) : DirichletTimeDependent() {}


        };

        class _DirichletTimeDependent {
public:

            /** Set kernels for constraint.
             *
             * @param[out] constraint Constraint for boundary condition.
             * @param[in] bc Dirichlet time-dependent boundary condition.
             * @param[in] solution Solution field.
             */
            static
            void setKernelConstraint(pylith::feassemble::ConstraintSpatialDB* constraint,
                                     const pylith::bc::DirichletTimeDependent& bc,
                                     const pylith::topology::Field& solution);

            static const char* pyreComponent;

        }; // _DirichletTimeDependent
        const char* _DirichletTimeDependent::pyreComponent = "dirichlettimedependent";

    } // bc
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletTimeDependent::DirichletTimeDependent(void) :
    _useInitial(true),
    _useRate(false),
    _useTimeHistory(false) {
    PyreComponent::setName(_DirichletTimeDependent::pyreComponent);
} // constructor


// ------------------------------------------------------------------------------------------------
// Factory for std::shared_ptr.
std::shared_ptr<pylith::bc::DirichletTimeDependent>
pylith::bc::DirichletTimeDependent::create(void) {
    return std::make_shared<DirichletTimeDependentWrap>();
}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletTimeDependent::~DirichletTimeDependent(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletTimeDependent::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    BoundaryCondition::deallocate();

    _dbTimeHistory.reset();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set indices of constrained degrees of freedom at each location.
void
pylith::bc::DirichletTimeDependent::setConstrainedDOF(const pylith::integer_array& dof) {
    PYLITH_COMPONENT_DEBUG("setConstrainedDOF(#dof="<<dof.size()<<")");

    _constrainedDOF = dof;
} // setConstrainedDOF


// ------------------------------------------------------------------------------------------------
// Get indices of constrained degrees of freedom.
const pylith::integer_array&
pylith::bc::DirichletTimeDependent::getConstrainedDOF(void) const {
    return _constrainedDOF;
} // getConstrainedDOF


// ------------------------------------------------------------------------------------------------
// Set time history database.
void
pylith::bc::DirichletTimeDependent::setTimeHistoryDB(const std::shared_ptr<spatialdata::spatialdb::TimeHistory>& th) {
    PYLITH_COMPONENT_DEBUG("setTimeHistoryDB(th"<<th<<")");

    _dbTimeHistory = th;
} // setTimeHistoryDB


// ------------------------------------------------------------------------------------------------
// Get time history database.
const spatialdata::spatialdb::TimeHistory&
pylith::bc::DirichletTimeDependent::getTimeHistoryDB(void) {
    return *_dbTimeHistory;
} // getTimeHistoryDB


// ------------------------------------------------------------------------------------------------
// Use initial value term in time history expression.
void
pylith::bc::DirichletTimeDependent::useInitial(const bool value) {
    PYLITH_COMPONENT_DEBUG("useInitial(value="<<value<<")");

    _useInitial = value;
} // useInitial


// ------------------------------------------------------------------------------------------------
// Get flag associated with using initial value term in time history expression.
bool
pylith::bc::DirichletTimeDependent::useInitial(void) const {
    return _useInitial;
} // useInitial


// ------------------------------------------------------------------------------------------------
// Use rate value term in time history expression.
void
pylith::bc::DirichletTimeDependent::useRate(const bool value) {
    PYLITH_COMPONENT_DEBUG("useRate(value="<<value<<")");

    _useRate = value;
} // useRate


// ------------------------------------------------------------------------------------------------
// Get flag associated with using rate value term in time history expression.
bool
pylith::bc::DirichletTimeDependent::useRate(void) const {
    return _useRate;
} // useRate


// ------------------------------------------------------------------------------------------------
// Use time history term in time history expression.
void
pylith::bc::DirichletTimeDependent::useTimeHistory(const bool value) {
    PYLITH_COMPONENT_DEBUG("useTimeHistory(value="<<value<<")");

    _useTimeHistory = value;
} // useTimeHistory


// ------------------------------------------------------------------------------------------------
// Get flag associated with using time history term in time history expression.
bool
pylith::bc::DirichletTimeDependent::useTimeHistory(void) const {
    return _useTimeHistory;
} // useTimeHistory


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::DirichletTimeDependent::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getName()<<")");

    if (!solution.hasSubfield(_subfieldName.c_str())) {
        std::ostringstream msg;
        msg << "Cannot constrain field '"<< _subfieldName
            << "' in component '" << PyreComponent::getIdentifier() << "'"
            << "; field is not in solution.";
        throw std::runtime_error(msg.str());
    } // if

    const topology::Field::SubfieldInfo& info = solution.getSubfieldInfo(_subfieldName.c_str());
    const int numComponents = info.description.numComponents;
    const int numConstrained = _constrainedDOF.size();
    for (int iConstrained = 0; iConstrained < numConstrained; ++iConstrained) {
        if (_constrainedDOF[iConstrained] >= numComponents) {
            std::ostringstream msg;
            msg << "Cannot constrain degree of freedom '" << _constrainedDOF[iConstrained] << "'"
                << " in component '" << PyreComponent::getIdentifier() << "'"
                << "; solution field '" << _subfieldName << "' contains only " << numComponents << " components.";
            throw std::runtime_error(msg.str());
        } // if
    } // for

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
std::shared_ptr<pylith::feassemble::Integrator>
pylith::bc::DirichletTimeDependent::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getName()<<") empty method");

    return nullptr;
} // createIntegrator


// ------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
std::vector<std::shared_ptr<pylith::feassemble::Constraint> >
pylith::bc::DirichletTimeDependent::createConstraints(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraints(solution="<<solution.getName()<<")");

    std::vector<std::shared_ptr<pylith::feassemble::Constraint> > constraintArray;
    std::shared_ptr<pylith::feassemble::ConstraintSpatialDB> constraint = pylith::feassemble::ConstraintSpatialDB::create(shared_from_this());assert(constraint);

    constraint->setSubfieldName(_subfieldName.c_str());
    constraint->setLabelName(getLabelName());
    constraint->setLabelValue(getLabelValue());
    constraint->setConstrainedDOF(_constrainedDOF);

    _DirichletTimeDependent::setKernelConstraint(constraint.get(), *this, solution);
    BoundaryCondition::_setKernelsDiagnosticField(constraint.get(), solution);

    constraintArray.resize(1);
    constraintArray[0] = constraint;

    PYLITH_METHOD_RETURN(constraintArray);
} // createConstraints


// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
std::shared_ptr<pylith::topology::Field>
pylith::bc::DirichletTimeDependent::createAuxiliaryField(const pylith::topology::Field& solution,
                                                         const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getName()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    std::shared_ptr<pylith::topology::Field> auxiliaryField = std::make_shared<pylith::topology::Field>(domainMesh);assert(auxiliaryField);
    auxiliaryField->setName("auxiliary field");

    assert(_subfieldFactory);
    assert(_normalizer);

    // :ATTENTION: The order of the factory methods must match the order of the auxiliary subfields in the FE kernels.

    _subfieldFactory->open(auxiliaryField, _normalizer);
    std::shared_ptr<pylith::topology::Field::Description> refDescription = std::make_shared<pylith::topology::Field::Description>(solution.getSubfieldInfo(_subfieldName.c_str()).description);
    _subfieldFactory->setRefDescription(refDescription);
    if (_useInitial) {
        _subfieldFactory->addSubfield(_subfieldFactory->initial_amplitude);
    } // if
    if (_useRate) {
        _subfieldFactory->addSubfield(_subfieldFactory->rate_amplitude);
        _subfieldFactory->addSubfield(_subfieldFactory->rate_start_time);
    } // _useRate
    if (_useTimeHistory) {
        _subfieldFactory->addSubfield(_subfieldFactory->time_history_amplitude);
        _subfieldFactory->addSubfield(_subfieldFactory->time_history_start_time);
        _subfieldFactory->addSubfield(_subfieldFactory->time_history_value);
        if (_dbTimeHistory) {
            _dbTimeHistory->open();
        } // if
    } // _useTimeHistory
    _subfieldFactory->close();
    refDescription.reset();

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    pylith::topology::FieldQuery fieldQuery(auxiliaryField);
    if (_useInitial) {
        _subfieldFactory->addSubfield(_subfieldFactory->initial_amplitude);
    } // if
    if (_useRate) {
        _subfieldFactory->addSubfield(_subfieldFactory->rate_amplitude);
        _subfieldFactory->addSubfield(_subfieldFactory->rate_start_time);
    } // _useRate
    if (_useTimeHistory) {
        _subfieldFactory->addSubfield(_subfieldFactory->time_history_amplitude);
        _subfieldFactory->addSubfield(_subfieldFactory->time_history_start_time);
    } // _useTimeHistory
    fieldQuery.open(_auxiliaryFieldDB, _normalizer->getLengthScale());
    fieldQuery.query();
    fieldQuery.close();

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying auxiliary field");
        auxiliaryField->view("Dirichlet auxiliary field");
    } // if

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::bc::DirichletTimeDependent::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                         const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateAuxiliaryField(auxiliaryField="<<auxiliaryField<<", t="<<t<<")");

    if (_useTimeHistory) {
        assert(_normalizer);
        const pylith::real timeScale = _normalizer->getTimeScale();
        TimeDependentOps::updateAuxiliaryField(auxiliaryField, t, timeScale, _dbTimeHistory);
    } // if

    PYLITH_METHOD_END;
} // updateAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Set kernels for computing constraint value.
void
pylith::bc::_DirichletTimeDependent::setKernelConstraint(pylith::feassemble::ConstraintSpatialDB* constraint,
                                                         const pylith::bc::DirichletTimeDependent& bc,
                                                         const topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_DirichletTimeDependent::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "setKernelConstraint(constraint="<<constraint<<", bc="<<typeid(bc).name()<<", solution="<<solution.getName()
          <<")" << pythia::journal::endl;

    PetscBdPointFunc bcKernel = nullptr;

    const pylith::topology::Field::VectorFieldEnum fieldType = solution.getSubfieldInfo(bc.getSubfieldName()).description.vectorFieldType;
    const bool isScalarField = fieldType == pylith::topology::Field::SCALAR;

    const int bitInitial = bc.useInitial() ? 0x1 : 0x0;
    const int bitRate = bc.useRate() ? 0x2 : 0x0;
    const int bitTimeHistory = bc.useTimeHistory() ? 0x4 : 0x0;
    const int bitUse = bitInitial | bitRate | bitTimeHistory;
    switch (bitUse) {
    case 0x1:
        bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::initial_scalar_boundary :
                   pylith::fekernels::TimeDependentFn::initial_vector_boundary;
        break;
    case 0x2:
        bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::rate_scalar_boundary :
                   pylith::fekernels::TimeDependentFn::rate_vector_boundary;
        break;
    case 0x4:
        bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::timeHistory_scalar_boundary :
                   pylith::fekernels::TimeDependentFn::timeHistory_vector_boundary;
        break;
    case 0x3:
        bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::initialRate_scalar_boundary :
                   pylith::fekernels::TimeDependentFn::initialRate_vector_boundary;
        break;
    case 0x5:
        bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::initialTimeHistory_scalar_boundary :
                   pylith::fekernels::TimeDependentFn::initialTimeHistory_vector_boundary;
        break;
    case 0x6:
        bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::rateTimeHistory_scalar_boundary :
                   pylith::fekernels::TimeDependentFn::rateTimeHistory_vector_boundary;
        break;
    case 0x7:
        bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::initialRateTimeHistory_scalar_boundary :
                   pylith::fekernels::TimeDependentFn::initialRateTimeHistory_vector_boundary;
        break;
    case 0x0: {
        pythia::journal::warning_t warning(_DirichletTimeDependent::pyreComponent);
        warning << pythia::journal::at(__HERE__)
                << "Dirichlet BC provides no constraints."
                << pythia::journal::endl;
        break;
    } // case 0x0
    default: {
        pythia::journal::error_t error(_DirichletTimeDependent::pyreComponent);
        error << pythia::journal::at(__HERE__)
              << "Unknown combination of flags for Dirichlet BC terms (useInitial="<<bc.useInitial()
              << ", useRate="<<bc.useRate()<<", useTimeHistory="<<bc.useTimeHistory()<<")."
              << pythia::journal::endl;
        throw std::logic_error("Unknown combination of flags for Dirichlet BC terms.");
    } // default
    } // switch

    assert(constraint);
    constraint->setKernelConstraint(bcKernel);

    PYLITH_METHOD_END;
} // setKernelConstraint


// End of file
