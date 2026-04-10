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

#include "pylith/fekernels/common/kernel.hh"
#include "pylith/fekernels/common/Matrix.hh"

#include <cassert>
#include <cstddef>


namespace pylith::fekernels::momentum {
    /// strain = 1/2 (F^T F − I) with F = I + ∇
    struct FiniteStrain {
        PYLITH_KERNEL static void compute(pylith::fekernels::Matrix3D& strain,
                                          const pylith::fekernels::Matrix3D& grad_u) {
            strain.zero();
            const size_t dim = strain.getDim();

            // F = I + grad_u
            pylith::fekernels::Matrix3D F;
            for (size_t i = 0; i < dim; i++) {
                for (size_t j = 0; j < dim; j++) {
                    F(i,j) = grad_u(i,j) + (i == j ? 1.0 : 0.0);
                } // for
            } // for

            // C = F^T F
            pylith::fekernels::Matrix3D C;
            C.zero();
            for (size_t i = 0; i < dim; i++) {
                for (size_t j = 0; j < dim; j++) {
                    for (size_t k = 0; k < dim; k++) {
                        C(i,j) += F(k,i) * F(k,j);
                    } // for
                } // for
            } // for

            // strain = 1/2 (C - I)
            for (size_t i = 0; i < dim; i++) {
                for (size_t j = 0; j < dim; j++) {
                    strain(i,j) = 0.5 * (C(i,j) - (i == j ? 1.0 : 0.0));
                } // for
            } // for
        } // compute

    }; // FiniteStrain

} // pylith::fekernels::momentum
