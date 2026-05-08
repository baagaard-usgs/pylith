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
#include "pylith/fekernels/common/Tensor2.hh"

#include <cstddef>


namespace pylith::fekernels::momentum {
    template<typename Dim, class SolutionLayout> class InfinitesimalStrain;
} // namespace


/// Infinitesimal strain
/// ε_ij = 1/2 (∂u_i/∂x_j + ∂u_j/∂x_i)
template<typename Dim, class SolutionLayout>
class pylith::fekernels::momentum::InfinitesimalStrain {
public:

    using SolutionUnpacked = typename SolutionLayout::template Unpacked<Dim>;

    PYLITH_KERNEL static void compute(pylith::fekernels::common::Tensor2<Dim>& strain,
                                      const SolutionUnpacked& solution);

}; // InfinitesimalStrain


#include "InfinitesimalStrain.icc"
