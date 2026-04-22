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

#include "pylith/fekernels/common/Kernel.hh"
#include "MomentumLayout.hh"

#include <cstddef>


namespace pylith::fekernels::momentum {
    template<typename Dim, class AuxiliaryLayout>
    class BodyForce;
} // namespce


/// f0 = f_i + ρ g_i
template<typename Dim, class AuxiliaryLayout>
class pylith::fekernels::momentum::BodyForce {
public:

    PYLITH_KERNEL static void f0(const pylith::integer cellDim,
                                 [[maybe_unused]] const pylith::integer numS,
                                 [[maybe_unused]] const pylith::integer numA,
                                 [[maybe_unused]] const pylith::integer sOff[],
                                 [[maybe_unused]] const pylith::integer sOff_x[],
                                 [[maybe_unused]] const pylith::scalar s[],
                                 [[maybe_unused]] const pylith::scalar s_t[],
                                 [[maybe_unused]] const pylith::scalar s_x[],
                                 const pylith::integer aOff[],
                                 const pylith::integer aOff_x[],
                                 const pylith::scalar a[],
                                 const pylith::scalar a_t[],
                                 const pylith::scalar a_x[],
                                 [[maybe_unused]] const pylith::real t,
                                 [[maybe_unused]] const pylith::real x[],
                                 [[maybe_unused]] const pylith::integer numConstants,
                                 [[maybe_unused]] const pylith::scalar constants[],
                                 pylith::scalar f0[]) noexcept;

}; // BodyForce


#include "BodyForce.icc"
