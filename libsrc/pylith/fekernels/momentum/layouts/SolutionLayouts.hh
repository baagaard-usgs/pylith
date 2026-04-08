// momentum/layouts/SolutionLayouts.hh
#pragma once
#include "../../common/Layout.hh"

namespace pylith::kernels::momentum::layout {
    /// Solution subfields for a purely mechanical problem.
    /// Subfield order: [displacement].
    struct MechanicalSolution {
        static constexpr int iDisp = 0;   ///< Displacement u_i
        static constexpr int size = 1;
    };

    /// Solution subfields for a hydromechanical (poromechanics) problem.
    /// Subfield order: [displacement, pore pressure].
    struct HydromechanicalSolution {
        static constexpr int iDisp = 0;  ///< Displacement u_i
        static constexpr int iPressure = 1;   ///< Pore pressure p
        static constexpr int size = 2;
    };

} // pylith::kernels::momentum::layout
