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
    template<typename Dim, class SolutionLayout> class DeviatoricStrain;
} // namespace


/// ε_dev = ε − 1/3 (tr ε) I
template<typename Dim, class SolutionLayout>
class pylith::fekernels::momentum::DeviatoricStrain {
public:

    /// Compute deviatoric strain given total strain.
    PYLITH_KERNEL static void compute(pylith::fekernels::common::Tensor2<Dim>& deviatoricStrain,
                                      const pylith::fekernels::common::Tensor2<Dim>& strain);

}; // DeviatoricStrain


#include "DeviatoricStrain.icc"
