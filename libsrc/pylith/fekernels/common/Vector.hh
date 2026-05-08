// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Vector storage and access.
 *
 * Provides local storage for a vector and access methods.
 */

#pragma once

#include "pylith/utils/types.hh"
#include "pylith/fekernels/common/Kernel.hh"

#include <cstddef>

namespace pylith::fekernels::common {
    template<typename Dim> class Vector;
} // namespace


/// Vector storage and access.
template <typename Dim>
class pylith::fekernels::common::Vector {
public:

    pylith::scalar _vector[Dim::value];

    PYLITH_KERNEL pylith::scalar& operator()(size_t i) noexcept;

    PYLITH_KERNEL pylith::scalar operator()(size_t i) const noexcept;

    PYLITH_KERNEL void zero(void) noexcept;

}; // Vector

#include "Vector.icc"
