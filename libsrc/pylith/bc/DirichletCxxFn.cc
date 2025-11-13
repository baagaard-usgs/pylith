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

#include "pylith/bc/DirichletCxxFn.hh" // implementation of object methods

#include "pylith/feassemble/ConstraintCxxFn.hh" // USES ConstraintBoundary
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
namespace pylith {
    namespace bc {
        class DirichletCxxFnWrap : public DirichletCxxFn {
public:

            DirichletCxxFnWrap(void) : DirichletCxxFn() {}


        };

    } // bc
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletCxxFn::DirichletCxxFn(void) :
    _fn(nullptr),
    _fnDot(nullptr) {
    PyreComponent::setName("dirichletuserfn");
} // constructor


// ------------------------------------------------------------------------------------------------
// Factory for std::shared_ptr.
std::shared_ptr<pylith::bc::DirichletCxxFn>
pylith::bc::DirichletCxxFn::create(void) {
    return std::make_shared<DirichletCxxFnWrap>();
}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletCxxFn::~DirichletCxxFn(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletCxxFn::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    BoundaryCondition::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set indices of constrained degrees of freedom at each location.
void
pylith::bc::DirichletCxxFn::setConstrainedDOF(const pylith::integer_array& dof) {
    PYLITH_COMPONENT_DEBUG("setConstrainedDOF(#dof="<<dof.size()<<")");

    _constrainedDOF = dof;
} // setConstrainedDOF


// ------------------------------------------------------------------------------------------------
// Get indices of constrained degrees of freedom.
const pylith::integer_array&
pylith::bc::DirichletCxxFn::getConstrainedDOF(void) const {
    return _constrainedDOF;
} // getConstrainedDOF


// ------------------------------------------------------------------------------------------------
// Set user function specifying field on boundary.
void
pylith::bc::DirichletCxxFn::setUserFn(PetscUserFieldFunc fn) {
    PYLITH_COMPONENT_DEBUG("setUserFn(fn="<<fn<<")");

    _fn = fn;
} // setUserFn


// ------------------------------------------------------------------------------------------------
// Get user function specifying field on boundary.
PetscUserFieldFunc
pylith::bc::DirichletCxxFn::getUserFn(void) const {
    return _fn;
} // getUserFn


// ------------------------------------------------------------------------------------------------
// Set user function specifying time derivative of field on boundary.
void
pylith::bc::DirichletCxxFn::setUserFnDot(PetscUserFieldFunc fn) {
    PYLITH_COMPONENT_DEBUG("setUserFnDot(fn="<<fn<<")");

    _fnDot = fn;
} // setUserFnDot


// ------------------------------------------------------------------------------------------------
// Get user function specifying time derivative of field on boundary.
PetscUserFieldFunc
pylith::bc::DirichletCxxFn::getUserFnDot(void) const {
    return _fnDot;
} // getUserFnDot


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::DirichletCxxFn::verifyConfiguration(const pylith::topology::Field& solution) const {
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
pylith::bc::DirichletCxxFn::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getName()<<") empty method");

    return nullptr;
} // createIntegrator


// ------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
std::vector<std::shared_ptr<pylith::feassemble::Constraint> >
pylith::bc::DirichletCxxFn::createConstraints(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraints(solution="<<solution.getName()<<")");

    std::vector<std::shared_ptr<pylith::feassemble::Constraint> > constraintArray;
    std::shared_ptr<pylith::feassemble::ConstraintCxxFn> constraint = pylith::feassemble::ConstraintCxxFn::create(shared_from_this());assert(constraint);
    constraint->setSubfieldName(_subfieldName.c_str());
    constraint->setLabelName(getLabelName());
    constraint->setLabelValue(getLabelValue());
    constraint->setConstrainedDOF(_constrainedDOF);
    constraint->setCxxFn(_fn);
    constraint->setCxxFnDot(_fnDot);

    constraintArray.resize(1);
    constraintArray[0] = constraint;
    PYLITH_METHOD_RETURN(constraintArray);
} // createConstraints


// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
std::shared_ptr<pylith::topology::Field>
pylith::bc::DirichletCxxFn::createAuxiliaryField(const pylith::topology::Field& solution,
                                                 const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getName()<<", domainMesh=)"<<typeid(domainMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(nullptr);
} // createAuxiliaryField


// End of file
