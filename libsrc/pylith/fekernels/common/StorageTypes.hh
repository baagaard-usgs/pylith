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
    template<typename Dim> struct Vector;
    template<typename Dim> struct Tensor2;
} // namespace


/// Vector storage and access.
template <typename Dim>
struct pylith::fekernels::common::Vector {
    pylith::scalar _vector[Dim::value];

    PYLITH_KERNEL pylith::scalar& operator()(size_t i) noexcept {
        assert(i < Dim::value);return _vector[i];
    } // operator()

    PYLITH_KERNEL pylith::scalar operator()(size_t i) const noexcept {
        assert(i < Dim::value);return _vector[i];
    } // operator()

    PYLITH_KERNEL void zero(void) noexcept {
        for (size_t i = 0; i < Dim::value; ++i) {
            _vector[i] = 0.0;
        } // for
    } // zero()

}; // Vector


/// Tensor2 storage and access.
template<typename Dim>
struct pylith::fekernels::common::Tensor2 {
    pylith::scalar _matrix[Dim::value][Dim::value];

    PYLITH_KERNEL pylith::scalar& operator()(size_t i,
                                             size_t j) noexcept {
        assert(i < Dim::value);assert(j < Dim::value);
        return _matrix[i][j];
    } // operator()

    PYLITH_KERNEL pylith::scalar operator()(size_t i,
                                            size_t j) const noexcept {
        assert(i < Dim::value);assert(j < Dim::value);
        return _matrix[i][j];
    } // operator()

    PYLITH_KERNEL void zero(void) noexcept {
        for (size_t i = 0; i < Dim::value; i++) {
            for (size_t j = 0; j < Dim::value; j++) {
                _matrix[i][j] = 0.0;
            } // for
        } // for
    } // zero

    PYLITH_KERNEL pylith::scalar trace(void) const noexcept {
        pylith::scalar value = 0.0;
        for (size_t i = 0; i < Dim::value; ++i) {
            value += _matrix[i][i];
        } // for
        return value;
    } // trace

}; // Tensor2
