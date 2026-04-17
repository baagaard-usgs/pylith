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


namespace pylith::fekernels::momentum {
    /// Flags indicating strain model
    enum class StrainFlags : uint32_t {
        INFINITESIMAL=0,
        FINITE=1,
    }; // StrainFlags

} // namespace
