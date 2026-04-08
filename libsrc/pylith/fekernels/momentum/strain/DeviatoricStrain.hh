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
#include "VolumetricStrain.hh"

#include <cassert>
#include <cstddef>


namespace pylith::fekernels::momentum {
    /// ε_dev = ε − 1/3 (tr ε) I
    struct DeviatoricStrain {
        PYLITH_KERNEL static void compute(pylith::fekernels::Matrix3D& deviatoricStrain,
                                          const pylith::fekernels::Matrix3D& strain) {
            const size_t dim = strain.getDim();

            const double volumetricStrain = VolumetricStrain::compute(strain);
            const double h = volumetricStrain / 3.0;

            for (size_t i = 0; i < dim; i++) {
                for (size_t j = 0; j < dim; j++) {
                    deviatoricStrain(i,j) = strain(i,j) - (i == j ? h : 0.0);
                } // for
            } // for
        } // compute

    }; // DeviatoricStrain

} // pylith::fekernels::momentum
