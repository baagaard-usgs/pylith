// common/Layout.hh
#pragma once

namespace pylith::layout {
    /// Wraps a Descriptor, adding Base to every index.
    /// Descriptor must provide: static constexpr int size, and one int per slot.
    /// Usage: inherit from BlockAt<MyDescriptor, startIndex>.
    template<class Descriptor, int Base>
    struct BlockAt : Descriptor {
        // Re-export shifted indices as a new type so kernels can rely on this.
        static constexpr int base = Base;

        // Validate at compile time that Base is non-negative.
        static_assert(Base >= 0, "Layout base index must be non-negative");
    };

    /// Utility: total size of two descriptors, for chaining bases.
    template<class D1, class D2>
    static constexpr int combined_size = D1::size + D2::size;

} // pylith::layout
