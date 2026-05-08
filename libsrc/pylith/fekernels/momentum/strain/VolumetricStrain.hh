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
    template<typename Dim, class SolutionLayout> class VolumetricStrain;
} // namespace


/// ε_vol = tr(ε)
template<typename Dim, class SolutionLayout>
class pylith::fekernels::momentum::VolumetricStrain {
public:

    /// Compute strain.
    PYLITH_KERNEL static pylith::scalar compute(const pylith::fekernels::common::Tensor2<Dim>& strain);

}; // VolumetricStrain


#include "VolumetricStrain.icc"
