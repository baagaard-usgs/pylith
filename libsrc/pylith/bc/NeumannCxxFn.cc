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

#include "pylith/bc/NeumannCxxFn.hh" // implementation of object methods

#include "pylith/fekernels/NeumannTimeDependent.hh" // USES NeumannTimeDepndent kernels

#include "pylith/feassemble/IntegratorBoundary.hh" // USES IntegratorBoundary
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/Mesh.hh" // USES Mesh

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
        class NeumannCxxFnWrap : public NeumannCxxFn {
public:

            NeumannCxxFnWrap(void) : NeumannCxxFn() {}


        };

        class _NeumannCxxFn {
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
                                    const pylith::bc::NeumannCxxFn& bc,
                                    const pylith::topology::Field& solution,
                                    const pylith::problems::Physics::FormulationEnum formulation);

            static const char* pyreComponent;

        }; // _NeumannCxxFn
        const char* _NeumannCxxFn::pyreComponent = "NeumannCxxFn";

    } // bc
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::NeumannCxxFn::NeumannCxxFn(void) :
    _fn(nullptr) {
    PyreComponent::setName("NeumannCxxFn");
} // constructor


// ------------------------------------------------------------------------------------------------
// Factory for std::shared_ptr.
std::shared_ptr<pylith::bc::NeumannCxxFn>
pylith::bc::NeumannCxxFn::create(void) {
    return std::make_shared<NeumannCxxFnWrap>();
}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::NeumannCxxFn::~NeumannCxxFn(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::NeumannCxxFn::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    BoundaryCondition::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set user function specifying field on boundary.
void
pylith::bc::NeumannCxxFn::setCxxFn(PetscBdPointFunc fn) {
    PYLITH_COMPONENT_DEBUG("setCxxFn(fn="<<fn<<")");

    _fn = fn;
} // setCxxFn


// ------------------------------------------------------------------------------------------------
// Get user function specifying field on boundary.
PetscBdPointFunc
pylith::bc::NeumannCxxFn::getCxxFn(void) const {
    return _fn;
} // getCxxFn


// ------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
std::shared_ptr<pylith::feassemble::Integrator>
pylith::bc::NeumannCxxFn::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getName()<<")");

    std::shared_ptr<pylith::feassemble::IntegratorBoundary> integrator = pylith::feassemble::IntegratorBoundary::create(shared_from_this());assert(integrator);
    integrator->setSubfieldName(getSubfieldName());
    integrator->setLabelName(getLabelName());
    integrator->setLabelValue(getLabelValue());

    _NeumannCxxFn::setKernelsResidual(integrator.get(), *this, solution, _formulation);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
std::vector<std::shared_ptr<pylith::feassemble::Constraint> >
pylith::bc::NeumannCxxFn::createConstraints(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraints(solution="<<solution.getName()<<") empty method");
    std::vector<std::shared_ptr<pylith::feassemble::Constraint> > constraintArray;

    PYLITH_METHOD_RETURN(constraintArray);
} // createConstraints


// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
std::shared_ptr<pylith::topology::Field>
pylith::bc::NeumannCxxFn::createAuxiliaryField(const pylith::topology::Field& solution,
                                               const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getName()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    std::shared_ptr<pylith::topology::Field> auxiliaryField = std::make_shared<pylith::topology::Field>(domainMesh);assert(auxiliaryField);
    auxiliaryField->setName("auxiliary field (not used)");

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying auxiliary field");
        auxiliaryField->view("Neumann auxiliary field");
    } // if

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::bc::_NeumannCxxFn::setKernelsResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                              const pylith::bc::NeumannCxxFn& bc,
                                              const topology::Field& solution,
                                              const pylith::problems::Physics::FormulationEnum formulation) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_NeumannCxxFn::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "setKernelsResidual(integrator="<<integrator<<", bc="<<typeid(bc).name()<<", solution="
          << solution.getName()<<")"
          << pythia::journal::endl;

    PetscBdPointFunc r0 = bc.getCxxFn();
    PetscBdPointFunc r1 = nullptr;

    std::vector<ResidualKernels> kernels(1);
    switch (formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        kernels[0] = ResidualKernels(bc.getSubfieldName(), pylith::feassemble::Integrator::LHS, r0, r1);
        break;
    case pylith::problems::Physics::DYNAMIC_IMEX:
    case pylith::problems::Physics::DYNAMIC:
        kernels[0] = ResidualKernels(bc.getSubfieldName(), pylith::feassemble::Integrator::RHS, r0, r1);
        break;
    default:
        PYLITH_JOURNAL_LOGICERROR("Unknown formulation for equations ("<<formulation<<").");
    } // switch

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // setKernelsResidual


// End of file
