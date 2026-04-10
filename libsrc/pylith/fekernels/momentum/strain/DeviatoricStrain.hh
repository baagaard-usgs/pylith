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
#include "VolumetricStrain.hh"

#include <cassert>
#include <cstddef>


namespace pylith::fekernels::momentum {
    template<int dim, class SolutionLayout> struct DeviatoricStrain;
} // namespace


/// ε_dev = ε − 1/3 (tr ε) I
template<int dim, class SolutionLayout>
struct pylith::fekernels::momentum::DeviatoricStrain {
    /// Compute deviatoric strain given total strain.
    PYLITH_KERNEL static void compute(pylith::fekernels::common::Matrix<dim>& deviatoricStrain,
                                      const pylith::fekernels::Matrix<dim>& strain) {
        const pylith::scalar volumetricStrain = VolumetricStrain<dim, SolutionLayout>::compute(strain);
        const pylith::scalar h = volumetricStrain / 3.0;

        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j < dim; j++) {
                deviatoricStrain(i,j) = strain(i,j) - (i == j ? h : 0.0);
            } // for
        } // for
    } // compute

}; // DeviatoricStrain
