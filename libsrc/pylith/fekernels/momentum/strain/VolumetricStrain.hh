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
    /// ε_vol = tr(ε)
    struct VolumetricStrain {
        PYLITH_KERNEL static pylith::scalar compute(const pylith::fekernels::Matrix3D& strain) {
            return strain.trace();
        } // compute

    }; // VolumetricStrain

} // pylith::fekernels::momentum
