// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** ScalarField provides access to a scalar field (value, time derivative, and gradient).
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
    class ScalarField;
}

class pylith::fekernels::common::ScalarField {
public:

    const pylith::scalar* value; ///< Value of field
    const pylith::scalar* value_t; ///< Time derivative of field
    const pylith::scalar* grad; ///< Gradient of field

    /// Get field value.
    PYLITH_KERNEL pylith::scalar operator()() const;

    /// Get time derivative of field.
    PYLITH_KERNEL pylith::scalar dot(void) const;

    /// Get gradient of field.
    PYLITH_KERNEL pylith::scalar gradient(size_t i) const;

}; // ScalarField

#include "ScalarField.icc"
