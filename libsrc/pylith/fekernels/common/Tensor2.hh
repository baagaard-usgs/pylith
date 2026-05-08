// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Tensor2 storage and access.
 *
 * Provides local storage for a rank 2 tensor (2D array) and access methods.
 */

#pragma once

#include "pylith/utils/types.hh"
#include "pylith/fekernels/common/Kernel.hh"

#include <cstddef>

namespace pylith::fekernels::common {
    template<typename Dim> class Tensor2;
} // namespace


/// Tensor2 storage and access.
template<typename Dim>
class pylith::fekernels::common::Tensor2 {
public:

    pylith::scalar _matrix[Dim::value][Dim::value];

    PYLITH_KERNEL pylith::scalar& operator()(size_t i,
                                             size_t j) noexcept;

    PYLITH_KERNEL pylith::scalar operator()(size_t i,
                                            size_t j) const noexcept;

    PYLITH_KERNEL void zero(void) noexcept;

    PYLITH_KERNEL pylith::scalar trace(void) const noexcept;

}; // Tensor2

#include "Tensor2.icc"
