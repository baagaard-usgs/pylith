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


#include <cstddef>


namespace pylith::fekernels::momentum {
    /// Flags for auxiliary subfields associated with momentum equation
    enum MomentumFlags : size_t {
        DEFAULT=0,
        BODY_FORCE=1 << 0,
        GRAVITY=1 << 1,
    }; // MomentumFlags

} // namespace
