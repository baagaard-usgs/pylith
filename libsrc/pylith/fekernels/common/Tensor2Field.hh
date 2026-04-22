// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Tensor2Field provides access to a rank-2 tensor field (value and time derivative).
 *
 * The object holds pointers to the storage already allocated.
 * The object is intended to be a member of an Unpacked object which provides named access
 * to the solution and auxiliary subfields.
 */

#pragma once

#include "pylith/utils/types.hh"
#include "kernel.hh"

#include <cassert>
#include <cstddef>

namespace pylith::fekernels::common {
    template<typename Dim> class Tensor2Field;
}

template <typename Dim>
class pylith::fekernels::common::Tensor2Field {
public:

    const pylith::scalar* value; ///< Value of field
    const pylith::scalar* value_t; ///< Time derivative of field
    const pylith::scalar* grad; ///< Gradient of field

    /// Get field value.
    PYLITH_KERNEL pylith::scalar operator()(size_t i,
                                            size_t j) const;

    /// Get time derivative of field.
    PYLITH_KERNEL pylith::scalar dot(size_t i,
                                     size_t j) const;

}; // Tensor2Field

#include "Tensor2Field.icc"
