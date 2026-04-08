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

#include "portability.hh"

#include <cassert>
#include <cstddef>

namespace pylith::fekernels {
    template<int dim>
    struct Matrix {
        double _matrix[dim][dim];
        PYLITH_KERNEL double& operator()(int i,
                                         int j) noexcept {
            return _matrix[i][j];
        } // operator()

        PYLITH_KERNEL double operator()(int i,
                                        int j) const noexcept {
            return _matrix[i][j];
        } // operator()

        PYLITH_KERNEL void zero(void) noexcept {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    _matrix[i][j] = 0.0;
                } // for
            } // for
        } // zero

        PYLITH_KERNEL size_t getDim(void) const noexcept {
            return dim;
        } // dim

        PYLITH_KERNEL double trace(void) const noexcept {
            double value = 0.0;
            for (size_t i = 0; i < dim; ++i) {
                value += _matrix[i][i];
            }
            return value;
        } // trace

    }; // Matrix
    using Matrix2D = Matrix<2>;
    using Matrix3D = Matrix<3>;

    PYLITH_KERNEL constexpr double
    delta(size_t i,
          size_t j) noexcept {
        return (i == j) ? 1.0 : 0.0;
    }


} // pylith::kernels
