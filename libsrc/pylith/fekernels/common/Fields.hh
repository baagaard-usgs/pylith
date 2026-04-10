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

#include "kernel.hh"

#include <cassert>
#include <cstddef>


namespace pylith::fekernels::common {
    struct ScalarField;
    template<size_t dim> struct VectorField;
}


struct pylith::fekernels::common::ScalarField {
    const pylith::scalar* value;
    const pylith::scalar* value_t;
    const pylith::scalar* grad;

    PYLITH_KERNEL pylith::scalar operator()() const {
        return value[0];
    }

    PYLITH_KERNEL pylith::scalar dot(void) const {
        return value_t[0];
    }

    PYLITH_KERNEL pylith::scalar grad_i(size_t i) const {
        return grad[i];
    }

}; // ScalarField


template <size_t dim>
struct pylith::fekernels::common::VectorField {
    const pylith::scalar* value;
    const pylith::scalar* value_t;
    const pylith::scalar* grad;

    PYLITH_KERNEL pylith::scalar operator()(size_t i) const {
        return value[i];
    } // operator()

    PYLITH_KERNEL pylith::scalar dot(size_t i) const {
        return value_t[i];
    } // dot()

    // grad[i*dim+j] = du_i/dx_j
    PYLITH_KERNEL pylith::scalar gradient(size_t i,
                                          size_t j) const {
        return grad[i*dim+j];
    } // grad_ij()

#if 0
    PYLITH_KERNEL pylith::scalar symmetricGradient(size_t i,
                                                   size_t j) const {
        return 0.5*(grad[i*dim+j] + grad[j*dim+i]);
    } // symmetricGradient()

    PYLITH_KERNEL pylith::scalar trace() const {
        pylith::scalar tr = 0.0;
        for (size_t i = 0; i < dim; i++) {
            tr += grad[i*dim+i];
        }
        return tr;
    } // trace()

    PYLITH_KERNEL pylith::scalar traceDot() const {
        pylith::scalar tr = 0.0;
        for (size_t i = 0; i < dim; i++) {
            tr += value_t[i*dim+i];
        }
        return tr;
    }

#endif
};
