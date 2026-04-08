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
    template<class StrainModel, class StressModel, class SolutionLayout, class AuxiliaryLayout>
    struct DivStress {
        PYLITH_KERNEL static void f1_domain(const PylithInt dim,
                                            const PylithInt numS,
                                            const PylithInt numA,
                                            const PylithInt sOff[],
                                            const PylithInt sOff_x[],
                                            const PylithScalar s[],
                                            const PylithScalar s_t[],
                                            const PylithScalar s_x[],
                                            const PylithInt aOff[],
                                            const PylithInt aOff_x[],
                                            const PylithScalar a[],
                                            const PylithScalar a_t[],
                                            const PylithScalar a_x[],
                                            const PylithReal t,
                                            const PylithReal x[],
                                            const PylithReal n[],
                                            const PylithInt numConstants,
                                            const PylithScalar constants[],
                                            PylithScalar f1[]) noexcept {
            pylith::fekernels::Matrix3D strain;
            pylith::fekernels::Matrix3D stress;
            pylith::fekernels::Matrix3D grad_u;
            strain.zero();
            stress.zero();
            stress.loadAuxiliary(a, aOff, AuxiliaryLayout);
            setGradU(grad_u, s_x[sOff_x[SolutionLayout::i_disp]], dim);
            grad_u.fromPointer(s_x[sOff_x[SolutionLayout::i_disp]]);

            StrainModel::compute(strain, grad_u);
            StressModel::compute(stress, strain);

            const size_t dim = stress.getDim();
            for (PylithInt i = 0; i < dim; ++i) {
                for (PylithInt j = 0; j < dim; ++j) {
                    f1[i*dim+j] = stress(i, j);
                } // for
            } // for

        } // f1_domain

    }; // DivStress

} // pylith::fekernels::momentum


/*
  Plane strain: ε_xz = ε_yz = 0, ε_zz = 0 (kinematic).
  We embed 2D grad_u (2x2) in 3D and compute σ using the 3D law:
    σ = 2 μ ε_dev + K tr(ε) I
  This yields a nonzero σ_zz from volumetric effects.
*/
struct ElasticIsotropic2DPlaneStrain {
  KERNEL_INLINE static void compute(const Mat2& grad_u_2d,
                                    const IsoElastic& mat,
                                    Mat3& sigma)
  {
    // Embed 2D grad_u into 3D
    Mat3 grad_u; grad_u.zero();
    grad_u(0,0) = grad_u_2d.v[0][0]; grad_u(0,1) = grad_u_2d.v[0][1];
    grad_u(1,0) = grad_u_2d.v[1][0]; grad_u(1,1) = grad_u_2d.v[1][1];
    // grad_u(*,2) and grad_u(2,*) remain 0 for plane strain

    // Reuse the 3D kernel
    ElasticIsotropic3D::compute(grad_u, mat, sigma);
  }
};

} // namespace fe