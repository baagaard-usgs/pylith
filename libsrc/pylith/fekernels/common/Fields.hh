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
    template<size_t dim> struct MatrixField;
}


struct pylith::fekernels::common::ScalarField {
    const pylith::scalar* value; ///< Value of field
    const pylith::scalar* value_t; ///< Time derivative of field
    const pylith::scalar* grad; ///< Gradient of field

    /// Get field value.
    PYLITH_KERNEL pylith::scalar operator()() const {
        return value[0];
    }

    /// Get time derivative of field.
    PYLITH_KERNEL pylith::scalar dot(void) const {
        return value_t[0];
    }

    /// Get gradient of field.
    PYLITH_KERNEL pylith::scalar gradient(size_t i) const {
        return grad[i];
    }

}; // ScalarField


template <size_t dim>
struct pylith::fekernels::common::VectorField {
    const pylith::scalar* value; ///< Value of field
    const pylith::scalar* value_t; ///< Time derivative of field
    const pylith::scalar* grad; ///< Gradient of field

    /// Get field value.
    PYLITH_KERNEL pylith::scalar operator()(size_t i) const {
        return value[i];
    } // operator()

    /// Get time derivative of field.
    PYLITH_KERNEL pylith::scalar dot(size_t i) const {
        return value_t[i];
    } // dot()

    /** Get gradient of field.
     *
     *  grad[i*dim+j] = du_i/dx_j
     */
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

template <size_t dim>
struct pylith::fekernels::common::MatrixField {
    const pylith::scalar* value; ///< Value of field
    const pylith::scalar* value_t; ///< Time derivative of field
    const pylith::scalar* grad; ///< Gradient of field

    /// Get field value.
    PYLITH_KERNEL pylith::scalar operator()(size_t i,
                                            size_t j) const {
        return value[i*dim+j];
    } // operator()

    /// Get time derivative of field.
    PYLITH_KERNEL pylith::scalar dot(size_t i,
                                     size_t j) const {
        return value_t[i*dim+j];
    } // dot()

#if 0
    /** Get gradient of field.
     *
     *  grad[i*dim+j] = du_i/dx_j
     */
    PYLITH_KERNEL pylith::scalar gradient(size_t i,
                                          size_t j) const {
        return grad[i*dim+j];
    } // grad_ij()

#endif

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
