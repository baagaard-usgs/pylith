// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** StorageTypes objects provide access to local storage of vectors, matrices, and tensors.
 *
 * The objects allocate local storage.
 */

#pragma once

#include "pylith/utils/types.hh"
#include "pylith/fekernels/common/kernel.hh"

#include <cstddef>

namespace pylith::fekernels::common {
    template<size_t dim> struct Vector;
    template<size_t dim> struct Matrix;
} // namespace


/// Vector storage and access.
template <size_t dim>
struct pylith::fekernels::common::Vector {
    pylith::scalar _vector[dim];

    PYLITH_KERNEL pylith::scalar& operator()(size_t i) noexcept {
        assert(i < dim);return _vector[i];
    } // operator()

    PYLITH_KERNEL pylith::scalar operator()(size_t i) const noexcept {
        assert(i < dim);return _vector[i];
    } // operator()

    PYLITH_KERNEL void zero(void) noexcept {
        for (size_t i = 0; i < dim; ++i) {
            _vector[i] = 0.0;
        } // for
    } // zero()

}; // Vector


/// Matrix storage and access.
template<size_t dim>
struct pylith::fekernels::common::Matrix {
    pylith::scalar _matrix[dim][dim];

    PYLITH_KERNEL pylith::scalar& operator()(size_t i,
                                             size_t j) noexcept {
        assert(i < dim);assert(j < dim);
        return _matrix[i][j];
    } // operator()

    PYLITH_KERNEL pylith::scalar operator()(size_t i,
                                            size_t j) const noexcept {
        assert(i < dim);assert(j < dim);
        return _matrix[i][j];
    } // operator()

    PYLITH_KERNEL void zero(void) noexcept {
        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j < dim; j++) {
                _matrix[i][j] = 0.0;
            } // for
        } // for
    } // zero

    PYLITH_KERNEL pylith::scalar trace(void) const noexcept {
        pylith::scalar value = 0.0;
        for (size_t i = 0; i < dim; ++i) {
            value += _matrix[i][i];
        } // for
        return value;
    } // trace

}; // Matrix
