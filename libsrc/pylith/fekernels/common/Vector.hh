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
    template <int dim>
    struct Vector {
        double _vector[dim];

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

        PYLITH_KERNEL size_t getDim(void) const noexcept {
            return dim;
        } // getDim

    };
    using Vector2D = Vector<2>;
    using Vector3D = Vector<3>;

} // pylith::fekernels
