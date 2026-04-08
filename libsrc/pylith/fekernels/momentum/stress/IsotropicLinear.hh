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

#include "pylith/fekernels/common/portability.hh"
#include "pylith/fekernels/common/Matrix.hh"

#include <cassert>
#include <cstddef>


namespace pylith::fekernels::momentum {
    template<class AuxiliaryLayout>
    struct IsotropicLinear {
        double _bulkModulus;
        double _shearModulus;

        // Local order of auxiliary subfields
        struct LocalLayout {
            static constexpr size_t i_bulkModulus = 0;
            static constexpr size_t i_shearModulus = 1;
            static constexpr size_t size = 2;
        };


        PYLITH_KERNEL void loadAuxiliary(const double* a,
                                         const int* aOff,
                                         pylith::fekernels::momentum::AuxiliaryLayout& layout) {
            _bulkModulus = a[aOff[layout::i_bulkModulus]];
            _shearModulus = a[aOff[layout::i_shearModulus]];
        } // loadAuxiliary

        PYLITH_KERNEL static void compute(pylith::fekernels::Matrix3D& stress,
                                          const pylith::fekernels::Matrix3D& strain) {
            const size_t dim = strain.getDim();

            meanStress(stress, _bulkModulus, strain);
            deviatoricStress(stress, _shearModulus, strain);

        } // compute

        /// C_ijkl = λ δ_ij δ_kl + μ(δ_ik δ_jl + δ_il δ_jk)
        PYLITH_KERNEL static void tangent(pylith::fekernels::Tensor43D& elasticityConstants) {
            const double lambda = _bulkModulus - 2.0/3.0 * _shearModulus;
            for (size_t i = 0; i < dim; ++i) {
                for (size_t j = 0; j < dim; ++j) {
                    for (size_t k = 0; k < dim; ++k) {
                        for (size_t l = 0; l < dim; ++l) {
                            elasticityConstants(i, j, k, l) =
                                lambda * delta(i, j) * delta(k, l) +
                                shearModulus * (delta(i, k) * delta(j, l) + (delta(i, l)* delta(j, k)));
                        } // for
                    } // for
                } // for
            } // for
        } // tangent

    }; // IsotropicLinear

} // pylith::fekernels::momentum
