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

#include "pylith/materials/RheologyPoroelasticity.hh" // implementation of object methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_DEBUG

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::RheologyPoroelasticity::RheologyPoroelasticity(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::RheologyPoroelasticity::~RheologyPoroelasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::RheologyPoroelasticity::deallocate(void) {
}


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::RheologyPoroelasticity::updateKernelConstants(pylith::real_array* kernelConstants,
                                                                 const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateKernelConstants(kernelConstants"<<kernelConstants<<", dt="<<dt<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables, implicit.
void
pylith::materials::RheologyPoroelasticity::addKernelsUpdateStateVarsImplicit(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                                                             const spatialdata::geocoords::CoordSys* coordsys,
                                                                             const bool _useStateVars) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVarsImplicit(kernels="<<kernels<<", coordsys="<<coordsys<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVarsImplicit


// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables, explicit.
void
pylith::materials::RheologyPoroelasticity::addKernelsUpdateStateVarsExplicit(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                                                             const spatialdata::geocoords::CoordSys* coordsys,
                                                                             const bool _useStateVars) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVarsExplicit(kernels="<<kernels<<", coordsys="<<coordsys<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVarsExplicit

//------------------------------------------------------------------------------------------------------------------------
// Unimplemented explicit kernel functions. Rheologies that use explicit time stepping should implement these.
PetscPointFunc 
pylith::materials::RheologyPoroelasticity::getKernelf0p_explicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0p_explicit is not implemented");

    return NULL;
} // getKernelf0p_explicit

PetscPointFunc 
pylith::materials::RheologyPoroelasticity::getKernelg0p(const spatialdata::geocoords::CoordSys* coordsys,
                                                        const bool _useBodyForce,
                                                        const bool _gravityField,
                                                        const bool _useSourceDensity) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg0p is not implemented");

    return NULL;
} // getKernelg0p

PetscPointFunc 
pylith::materials::RheologyPoroelasticity::getKernelg1p_explicit(const spatialdata::geocoords::CoordSys* coordsys,
                                                                 const bool _gravityField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1p_explicit is not implemented");

    return NULL;

} // getKernelg1p_explicit

PetscPointFunc 
pylith::materials::RheologyPoroelasticity::getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1v_explicit is not implemented");

    return NULL;
} // getKernelg1v_explicit



// End of file
