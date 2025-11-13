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

#include "pylith/bc/AbsorbingDampers.hh" // implementation of object methods

#include "pylith/bc/SubfieldFactory.hh" // USES SubfieldFactory

#include "pylith/fekernels/AbsorbingDampers.hh" // USES AbsorbingDampers kernels

#include "pylith/feassemble/IntegratorBoundary.hh" // USES IntegratorBoundary
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES typeid

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorBoundary::ResidualKernels ResidualKernels;

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace bc {
        class AbsorbingDampersWrap : public AbsorbingDampers {
public:

            AbsorbingDampersWrap(void) : AbsorbingDampers() {}


        };

        class _AbsorbingDampers {
public:

            /** Set kernels for RHS residual.
             *
             * @param[out] integrator Integrator for boundary condition.
             * @param[in] bc Absorbing dampers boundary condition.
             * @param[in] solution Solution field.
             */
            static
            void setKernelsResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                    const pylith::bc::AbsorbingDampers& bc,
                                    const pylith::topology::Field& solution);

            static const char* pyreComponent;
        };
        const char* _AbsorbingDampers::pyreComponent = "absorbingdampers";

    } // bc
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::AbsorbingDampers::AbsorbingDampers(void) {
    PyreComponent::setName(_AbsorbingDampers::pyreComponent);

    _subfieldName = "velocity";
} // constructor


// ------------------------------------------------------------------------------------------------
// Factory for std::shared_ptr.
std::shared_ptr<pylith::bc::AbsorbingDampers>
pylith::bc::AbsorbingDampers::create(void) {
    return std::make_shared<AbsorbingDampersWrap>();
}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::AbsorbingDampers::~AbsorbingDampers(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::AbsorbingDampers::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    BoundaryCondition::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::AbsorbingDampers::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getName()<<")");

    BoundaryCondition::verifyConfiguration(solution);

    const pylith::topology::Field::SubfieldInfo& info = solution.getSubfieldInfo(_subfieldName.c_str());
    if (pylith::topology::Field::VECTOR != info.description.vectorFieldType) {
        std::ostringstream msg;
        msg << "Absorbing boundary condition cannot be applied to non-vector field '"<< _subfieldName << "' in solution.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
std::shared_ptr<pylith::feassemble::Integrator>
pylith::bc::AbsorbingDampers::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getName()<<")");

    std::shared_ptr<pylith::feassemble::IntegratorBoundary> integrator = pylith::feassemble::IntegratorBoundary::create(shared_from_this());assert(integrator);
    integrator->setSubfieldName(getSubfieldName());
    integrator->setLabelName(getLabelName());
    integrator->setLabelValue(getLabelValue());

    _AbsorbingDampers::setKernelsResidual(integrator.get(), *this, solution);
    BoundaryCondition::_setKernelsDiagnosticField(integrator.get(), solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
std::vector<std::shared_ptr<pylith::feassemble::Constraint> >
pylith::bc::AbsorbingDampers::createConstraints(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraints(solution="<<solution.getName()<<") empty method");
    std::vector<std::shared_ptr<pylith::feassemble::Constraint> > constraints;

    PYLITH_METHOD_RETURN(constraints);
} // createConstraints


// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
std::shared_ptr<pylith::topology::Field>
pylith::bc::AbsorbingDampers::createAuxiliaryField(const pylith::topology::Field& solution,
                                                   const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getName()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    std::shared_ptr<pylith::topology::Field> auxiliaryField(new pylith::topology::Field(domainMesh));assert(auxiliaryField);
    auxiliaryField->setName("auxiliary field");

    assert(_subfieldFactory);
    assert(_normalizer);

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.
    _subfieldFactory->open(auxiliaryField, _normalizer);
    _subfieldFactory->addSubfield(_subfieldFactory->density); // 0
    _subfieldFactory->addSubfield(_subfieldFactory->vp); // 1
    _subfieldFactory->addSubfield(_subfieldFactory->vs); // 2
    _subfieldFactory->close();

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    pylith::topology::FieldQuery fieldQuery(auxiliaryField);
    fieldQuery.addSubfield(_subfieldFactory->density);
    fieldQuery.addSubfield(_subfieldFactory->vp);
    fieldQuery.addSubfield(_subfieldFactory->vs);
    fieldQuery.open(_auxiliaryFieldDB, _normalizer->getLengthScale());
    fieldQuery.query();
    fieldQuery.close();

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::bc::_AbsorbingDampers::setKernelsResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                                  const pylith::bc::AbsorbingDampers& bc,
                                                  const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_AbsorbingDampers::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "_AbsorbingDampers::_setKernelsRHSResidual(integrator="<<integrator<<", bc="<<typeid(bc).name()
          <<", solution="<<solution.getName()<<")"<<pythia::journal::endl;

    PetscBdPointFunc g0 = pylith::fekernels::AbsorbingDampers::g0;
    PetscBdPointFunc g1 = nullptr;

    std::vector<ResidualKernels> kernels(1);
    kernels[0] = ResidualKernels(bc.getSubfieldName(), pylith::feassemble::Integrator::RHS, g0, g1);

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsRHSResidual


// End of file
