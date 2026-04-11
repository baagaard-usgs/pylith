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
#include "pylith/fekernels/common/StorageTypes.hh"

#include <cassert>
#include <cstddef>

namespace pylith::fekernels::momentum {
    template<size_t dim, class SolutionLayout> struct VolumetricStrain;
} // namespace


/// ε_vol = tr(ε)
template<size_t dim, class SolutionLayout>
struct pylith::fekernels::momentum::VolumetricStrain {
    /// @brief  Compute strain.
    PYLITH_KERNEL static pylith::scalar compute(const pylith::fekernels::common::Matrix<dim>& strain) {
        return strain.trace();
    } // compute

}; // VolumetricStrain
