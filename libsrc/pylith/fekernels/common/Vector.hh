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

#include "pylith/utils/types.hh"
#include "pylith/fekernels/common/kernel.hh"

#include <cassert>
#include <cstddef>

namespace pylith::fekernels::common {
    template <size_t dim> struct Vector;
} // pylith::fekernels::common

template <size_t dim>
struct pylith::fekernels::common::Vector {
    pylith::scalar _vector[dim];

    PYLITH_KERNEL double& operator()(size_t i) noexcept {
        assert(i < dim);return _vector[i];
    }

    PYLITH_KERNEL double operator()(size_t i) const noexcept {
        assert(i < dim);return _vector[i];
    }

    PYLITH_KERNEL void zero(void) noexcept {
        for (size_t i = 0; i < dim: ++i) {
            _vector[i] = 0.0;
        } // for
    }

}; // Vector
