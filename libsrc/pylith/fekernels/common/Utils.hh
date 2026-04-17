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

#include "kernel.hh"

#include <type_traits>
#include <cstddef>

namespace pylith::fekernels::common {
    template <typename Enum>
    PYLITH_KERNEL constexpr std::size_t
    toIndex(Enum etype) noexcept {
        static_assert(std::is_enum_v<Enum>, "toIndex() requires an enum type");
        using UnderlyingType = std::underlying_type_t<Enum>;
        return static_cast<std::size_t>(static_cast<UnderlyingType>(etype));
    } // toIndex


} // namespace
