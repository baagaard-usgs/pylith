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


namespace pylith::fekernels::common {
    /// Optional subfields helper using empty base class optimization.
    /// When `present` is false, OptionalMember<false, T> has zero size.
    /// When `present` is true, OptionalMember<true, T> contains a T member.
    template <bool present, typename T>
    struct OptionalMember {};

    template <typename T>
    struct OptionalMember<true, T> { T member; };

    /// Helper to compute field count contribution (0 or 1 based on condition).
    constexpr int
    addIf(bool condition) noexcept {
        return condition ? 1 : 0;
    }


} // namespace pylith::fekernels::common
