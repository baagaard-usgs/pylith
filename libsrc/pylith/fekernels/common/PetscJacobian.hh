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

#include "pylith/fekernels/common/Kernel.hh"

#include <cstddef>


namespace pylith::fekernels::common {
    template<typename Dim> class PetscJacobian;
}

template<typename Dim>
class pylith::fekernels::common::PetscJacobian {
public:

    /** Get index in J3 Jacobian term.
     *
     * @param[in] f Trial function component
     * @param[in] df Derivative component of trial function
     * @param[in] g Basis function component
     * @param[in] dg Derivative component of basis function
     * @returns Index in J3 Jacobian term.
     */
    PYLITH_KERNEL static size_t index3(const size_t f,
                                       const size_t df,
                                       const size_t g,
                                       const size_t dg);

}; // Jacobian

#include "PetscJacobian.icc"
