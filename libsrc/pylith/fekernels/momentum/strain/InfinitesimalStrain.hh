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
    template<int dim, class SolutionLayout> struct InfinitesimalStrain;
} // namespace


/// Infinitesimal strain
/// ε_ij = 1/2 (∂u_i/∂x_j + ∂u_j/∂x_i)
template<int dim, class SolutionLayout>
struct pylith::fekernels::momentum::InfinitesimalStrain {
    PYLITH_KERNEL static void compute(pylith::fekernels::common::Matrix<dim>& strain,
                                      const SolutionLayout& solution) {
        strain.zero();
        const pylith::fekernels::common::VectorField<dim>& disp = solution.displacement;
        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j < dim; j++) {
                strain(i, j) = 0.5 * (disp.gradient(i, j) + disp.gradient(j, i));
            } // for
        } // for
    } // compute

}; // InfinitesimalStrain
