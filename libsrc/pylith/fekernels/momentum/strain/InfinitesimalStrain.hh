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

#include "pylith/fekernels/common/portability.hh"
#include "pylith/fekernels/common/Matrix.hh"

#include <cassert>
#include <cstddef>


namespace pylith::fekernels::momentum {
    /// ε_ij = 1/2 (∂u_i/∂x_j + ∂u_j/∂x_i)
    struct InfinitesimalStrain {
        PYLITH_KERNEL static void compute(pylith::fekernels::Matrix3D& strain,
                                          const pylith::fekernels::Matrix3D& grad_u) {
            strain.zero();
            const size_t dim = strain.getDim();
            for (size_t i = 0; i < dim; i++) {
                for (size_t j = 0; j < dim; j++) {
                    strain(i, j) = 0.5 * (grad_u(i, j) + grad_u(j, i));
                } // for
            } // for
        } // compute

    }; // InfinitesimalStrain

} // pylith::fekernels::momentum
