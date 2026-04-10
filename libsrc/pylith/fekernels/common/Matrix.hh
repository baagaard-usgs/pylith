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

#include <cstddef>

namespace pylith::fekernels::common {
    template<size_t dim> struct Matrix;
} // namespace


/// Matrix storage and access.
template<size_t dim>
struct pylith::fekernels::common::Matrix {
    pylith::scalar _matrix[dim][dim];

    PYLITH_KERNEL double& operator()(int i,
                                     int j) noexcept {
        assert(i < dim);assert(j < dim);
        return _matrix[i][j];
    } // operator()

    PYLITH_KERNEL double operator()(int i,
                                    int j) const noexcept {
        assert(i < dim);assert(j < dim);
        return _matrix[i][j];
    } // operator()

    PYLITH_KERNEL void zero(void) noexcept {
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                _matrix[i][j] = 0.0;
            } // for
        } // for
    } // zero

    PYLITH_KERNEL double trace(void) const noexcept {
        double value = 0.0;
        for (size_t i = 0; i < dim; ++i) {
            value += _matrix[i][i];
        } // for
        return value;
    } // trace

}; // Matrix
