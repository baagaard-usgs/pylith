// momentum/stress/IsotropicLinear.hh  (layout portion)
#pragma once
#include "../../common/Layout.hh"

namespace pylith::kernels::momentum {
    struct IsotropicLinear {
        /// Auxiliary subfields consumed by this rheology.
        /// 0-relative; composed layouts shift these via BlockAt<>.
        struct Desc {
            static constexpr int iK = 0;  ///< Bulk modulus K
            static constexpr int iMu = 1;  ///< Shear modulus μ
            static constexpr int size = 2;
        };

        // ... static compute / computeTangent methods (unchanged from before)
    };

} // pylith::kernels::momentum
