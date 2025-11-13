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

#include "pylith/bc/NeumannTimeDependent.hh" // implementation of object methods

#include "pylith/bc/SubfieldFactory.hh" // USES SubfieldFactory
#include "pylith/bc/TimeDependentOps.hh" // USES TimeDependentOps

#include "pylith/fekernels/NeumannTimeDependent.hh" // USES NeumannTimeDependent kernels

#include "pylith/feassemble/IntegratorBoundary.hh" // USES IntegratorBoundary
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorBoundary::ResidualKernels ResidualKernels;

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace bc {
        class NeumannTimeDependentWrap : public NeumannTimeDependent {
public:

            NeumannTimeDependentWrap(void) : NeumannTimeDependent() {}


        };

        class _NeumannTimeDependent {
public:

            /** Set kernels for RHS residual.
             *
             * @param[out] integrator Integrator for boundary condition.
             * @param[in] bc Neumann time-dependent boundary condition.
             * @param[in] solution Solution field.
             * @param[in] formulation Formulation for equations.
             */
            static
            void setKernelsResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                    const pylith::bc::NeumannTimeDependent& bc,
                                    const pylith::topology::Field& solution,
                                    const pylith::problems::Physics::FormulationEnum formulation);

            static const char* pyreComponent;

        }; // _NeumannTimeDependent
        const char* _NeumannTimeDependent::pyreComponent = "neumanntimedependent";

    } // bc
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::NeumannTimeDependent::NeumannTimeDependent(void) :
    _scaleName("pressure"),
    _useInitial(true),
    _useRate(false),
    _useTimeHistory(false) {
    PyreComponent::setName(_NeumannTimeDependent::pyreComponent);
} // constructor


// ------------------------------------------------------------------------------------------------
// Factory for std::shared_ptr.
std::shared_ptr<pylith::bc::NeumannTimeDependent>
pylith::bc::NeumannTimeDependent::create(void) {
    return std::make_shared<NeumannTimeDependentWrap>();
}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::NeumannTimeDependent::~NeumannTimeDependent(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::NeumannTimeDependent::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    BoundaryCondition::deallocate();

    _dbTimeHistory.reset();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set time history database.
void
pylith::bc::NeumannTimeDependent::setTimeHistoryDB(const std::shared_ptr<spatialdata::spatialdb::TimeHistory>& th) {
    PYLITH_COMPONENT_DEBUG("setTimeHistoryDB(th"<<th<<")");

    _dbTimeHistory = th;
} // setTimeHistoryDB


// ------------------------------------------------------------------------------------------------
// Get time history database.
const std::shared_ptr<spatialdata::spatialdb::TimeHistory>&
pylith::bc::NeumannTimeDependent::getTimeHistoryDB(void) {
    return _dbTimeHistory;
} // getTimeHistoryDB


// ------------------------------------------------------------------------------------------------
// Use initial value term in time history expression.
void
pylith::bc::NeumannTimeDependent::useInitial(const bool value) {
    PYLITH_COMPONENT_DEBUG("useInitial(value="<<value<<")");

    _useInitial = value;
} // useInitial


// ------------------------------------------------------------------------------------------------
// Get flag associated with using initial value term in time history expression.
bool
pylith::bc::NeumannTimeDependent::useInitial(void) const {
    return _useInitial;
} // useInitial


// ------------------------------------------------------------------------------------------------
// Use rate value term in time history expression.
void
pylith::bc::NeumannTimeDependent::useRate(const bool value) {
    PYLITH_COMPONENT_DEBUG("useRate(value="<<value<<")");

    _useRate = value;
} // useRate


// ------------------------------------------------------------------------------------------------
// Get flag associated with using rate value term in time history expression.
bool
pylith::bc::NeumannTimeDependent::useRate(void) const {
    return _useRate;
} // useRate


// ------------------------------------------------------------------------------------------------
// Use time history term in time history expression.
void
pylith::bc::NeumannTimeDependent::useTimeHistory(const bool value) {
    PYLITH_COMPONENT_DEBUG("useTimeHistory(value="<<value<<")");

    _useTimeHistory = value;
} // useTimeHistory


// ------------------------------------------------------------------------------------------------
// Get flag associated with using time history term in time history expression.
bool
pylith::bc::NeumannTimeDependent::useTimeHistory(void) const { // useTimeHistory
    return _useTimeHistory;
} // useTimeHistory


// ------------------------------------------------------------------------------------------------
// Name of scale associated with Neumann boundary condition (e.g., pressure for elasticity).
void
pylith::bc::NeumannTimeDependent::setScaleName(const char* value) {
    PYLITH_COMPONENT_DEBUG("setScaleName(value"<<value<<")");

    if (( value == std::string("length")) ||
        ( value == std::string("time")) ||
        ( value == std::string("pressure")) ||
        ( value == std::string("density")) ||
        ( value == std::string("pressure")) ) {
        _scaleName = value;
    } else {
        std::ostringstream msg;
        msg << "Unknown name of scale ("<<value<<") for Neumann boundary condition '" << getLabelName() << "'.";
        throw std::runtime_error(msg.str());
    } // if
} // setScaleName


// ------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
std::shared_ptr<pylith::feassemble::Integrator>
pylith::bc::NeumannTimeDependent::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getName()<<")");

    std::shared_ptr<pylith::feassemble::IntegratorBoundary> integrator = pylith::feassemble::IntegratorBoundary::create(shared_from_this());assert(integrator);
    integrator->setSubfieldName(getSubfieldName());
    integrator->setLabelName(getLabelName());
    integrator->setLabelValue(getLabelValue());

    _NeumannTimeDependent::setKernelsResidual(integrator.get(), *this, solution, _formulation);
    BoundaryCondition::_setKernelsDiagnosticField(integrator.get(), solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
std::vector<std::shared_ptr<pylith::feassemble::Constraint> >
pylith::bc::NeumannTimeDependent::createConstraints(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraints(solution="<<solution.getName()<<") empty method");
    std::vector<std::shared_ptr<pylith::feassemble::Constraint> > constraintArray;

    PYLITH_METHOD_RETURN(constraintArray);
} // createConstraints


// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
std::shared_ptr<pylith::topology::Field>
pylith::bc::NeumannTimeDependent::createAuxiliaryField(const pylith::topology::Field& solution,
                                                       const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getName()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    std::shared_ptr<pylith::topology::Field> auxiliaryField = std::make_shared<pylith::topology::Field>(domainMesh);assert(auxiliaryField);
    auxiliaryField->setName("auxiliary field");

    assert(_subfieldFactory);
    assert(_normalizer);
    std::shared_ptr<pylith::topology::Field::Description> refDescription = std::make_shared<pylith::topology::Field::Description>(solution.getSubfieldInfo(_subfieldName.c_str()).description);
    assert(refDescription);
    if (_scaleName == std::string("pressure")) {
        refDescription->scale = _normalizer->getPressureScale();
    } else if (_scaleName == std::string("velocity")) {
        refDescription->scale = sqrt(_normalizer->getPressureScale() / _normalizer->getDensityScale());
    } else if (_scaleName == std::string("length")) {
        refDescription->scale = _normalizer->getLengthScale();
    } else if (_scaleName == std::string("time")) {
        refDescription->scale = sqrt(_normalizer->getDensityScale() / _normalizer->getPressureScale()) * _normalizer->getLengthScale();
    } else if (_scaleName == std::string("density")) {
        refDescription->scale = _normalizer->getDensityScale();
    } else {
        std::ostringstream msg;
        msg << "Unknown name of scale ("<<_scaleName<<") for Neumann boundary condition for '" << getLabelName() << "'.";
        PYLITH_COMPONENT_ERROR(msg.str());
        throw std::logic_error(msg.str());
    } // if/else

    // :ATTENTION: The order of the factory methods must match the order of the auxiliary subfields in the FE kernels.
    _subfieldFactory->open(auxiliaryField, _normalizer);
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
        auxiliaryField->view("Neumann auxiliary field");
    } // if

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Update auxiliary subfields at beginning of time step.
void
pylith::bc::NeumannTimeDependent::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
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
// Set kernels for residual.
void
pylith::bc::_NeumannTimeDependent::setKernelsResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                                      const pylith::bc::NeumannTimeDependent& bc,
                                                      const topology::Field& solution,
                                                      const pylith::problems::Physics::FormulationEnum formulation) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_NeumannTimeDependent::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "setKernelsResidual(integrator="<<integrator<<", bc="<<typeid(bc).name()<<", solution="
          << solution.getName()<<")"
          << pythia::journal::endl;

    const pylith::topology::Field::VectorFieldEnum fieldType = solution.getSubfieldInfo(bc.getSubfieldName()).description.vectorFieldType;
    const bool isScalarField = fieldType == pylith::topology::Field::SCALAR;

    const int bitInitial = bc.useInitial() ? 0x1 : 0x0;
    const int bitRate = bc.useRate() ? 0x2 : 0x0;
    const int bitTimeHistory = bc.useTimeHistory() ? 0x4 : 0x0;
    const int bitUse = bitInitial | bitRate | bitTimeHistory;

    PetscBdPointFunc r0 = nullptr;
    PetscBdPointFunc r1 = nullptr;
    switch (bitUse) {
    case 0x1:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::f0_initial_scalar :
             pylith::fekernels::NeumannTimeDependent::f0_initial_vector;
        break;
    case 0x2:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::f0_rate_scalar :
             pylith::fekernels::NeumannTimeDependent::f0_rate_vector;
        break;
    case 0x4:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::f0_timeHistory_scalar :
             pylith::fekernels::NeumannTimeDependent::f0_timeHistory_vector;
        break;
    case 0x3:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::f0_initialRate_scalar :
             pylith::fekernels::NeumannTimeDependent::f0_initialRate_vector;
        break;
    case 0x5:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::f0_initialTimeHistory_scalar :
             pylith::fekernels::NeumannTimeDependent::f0_initialTimeHistory_vector;
        break;
    case 0x6:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::f0_rateTimeHistory_scalar :
             pylith::fekernels::NeumannTimeDependent::f0_rateTimeHistory_vector;
        break;
    case 0x7:
        r0 = (isScalarField) ? pylith::fekernels::NeumannTimeDependent::f0_initialRateTimeHistory_scalar :
             pylith::fekernels::NeumannTimeDependent::f0_initialRateTimeHistory_vector;
        break;
    case 0x0: {
        pythia::journal::warning_t warning(_NeumannTimeDependent::pyreComponent);
        warning << pythia::journal::at(__HERE__)
                << "Neumann time-dependent BC provides no values."
                << pythia::journal::endl;
        break;
    } // case 0x0
    default: {
        PYLITH_JOURNAL_LOGICERROR("Unknown combination of flags for Neumann BC terms (useInitial="
                                  <<bc.useInitial()
                                  << ", useRate="<<bc.useRate()<<", useTimeHistory="<<bc.useTimeHistory()<<").");
    } // default
    } // switch

    std::vector<ResidualKernels> kernels(1);
    switch (formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        kernels[0] = ResidualKernels(bc.getSubfieldName(), pylith::feassemble::Integrator::LHS, r0, r1);
        break;
    case pylith::problems::Physics::DYNAMIC_IMEX:
    case pylith::problems::Physics::DYNAMIC:
        kernels[0] = ResidualKernels(bc.getSubfieldName(), pylith::feassemble::Integrator::RHS, r0, r1);
        break;
    default: {
        PYLITH_JOURNAL_LOGICERROR("Unknown formulation for equations ("<<formulation<<").");
    } // default
    } // switch

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // setKernelsResidual


// End of file
