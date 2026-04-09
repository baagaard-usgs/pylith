#include <petsc.h>

// ============================================================
// Field accessor types
// ============================================================

template <PylithInt dim>
struct VectorField {
    const PylithScalar* value;
    const PylithScalar* value_t;
    const PylithScalar* grad;

    PETSC_DEVICE inline PylithScalar
    operator()(PylithInt i) const { return value[i]; }

    PETSC_DEVICE inline PylithScalar
    dot(PylithInt i) const { return value_t[i]; }

    // grad[i*dim+j] = du_i/dx_j
    PETSC_DEVICE inline PylithScalar
    grad_ij(PylithInt i, PylithInt j) const { return grad[i*dim+j]; }

    PETSC_DEVICE inline PylithScalar
    symGrad(PylithInt i, PylithInt j) const {
        return 0.5*(grad[i*dim+j] + grad[j*dim+i]);
    }

    PETSC_DEVICE inline PylithScalar
    trace() const {
        PylithScalar tr = 0.0;
        for (PylithInt i = 0; i < dim; i++) tr += grad[i*dim+i];
        return tr;
    }

    PETSC_DEVICE inline PylithScalar
    traceDot() const {
        PylithScalar tr = 0.0;
        for (PylithInt i = 0; i < dim; i++) tr += value_t[i*dim+i];
        return tr;
    }
};

struct ScalarField {
    const PylithScalar* value;
    const PylithScalar* value_t;
    const PylithScalar* grad;

    PETSC_DEVICE inline PylithScalar
    operator()() const { return value[0]; }

    PETSC_DEVICE inline PylithScalar
    dot() const { return value_t[0]; }

    PETSC_DEVICE inline PylithScalar
    grad_i(PylithInt i) const { return grad[i]; }
};

// ============================================================
// IsotropicLinearElasticity
// ============================================================

struct IsotropicLinearElasticity {

    // ----------------------------------------------------------
    // Solution layout
    // ----------------------------------------------------------

    struct Solution {
        enum Fields : PylithInt { DISPLACEMENT = 0, NUM_FIELDS };

        template <PylithInt dim>
        struct Unpacked {
            VectorField<dim> displacement;
        };

        template <PylithInt dim>
        PETSC_DEVICE static inline Unpacked<dim>
        unpack(const PylithInt sOff[],
               const PylithInt sOff_x[],
               const PylithScalar s[],
               const PylithScalar s_t[],
               const PylithScalar s_x[]) {
            return { { &s  [sOff  [DISPLACEMENT]],
                       &s_t[sOff  [DISPLACEMENT]],
                       &s_x[sOff_x[DISPLACEMENT]] } };
        }
    };

    // ----------------------------------------------------------
    // Auxiliary layout
    //
    // Field order must match the order fields are added to the
    // PyLith auxiliary factory:
    //   0: density    (scalar)
    //   1: vs         (scalar) -- stored, converted to mu on unpack
    //   2: vp         (scalar) -- stored, converted to lambda on unpack
    // ----------------------------------------------------------

    struct Auxiliary {
        enum Fields : PylithInt {
            DENSITY = 0,
            VS      = 1,
            VP      = 2,
            NUM_FIELDS
        };

        struct Unpacked {
            // Stored as named scalars — no local tensor, pure registers
            PylithScalar density;
            PylithScalar mu;        // = density * vs^2
            PylithScalar lambda;    // = density * (vp^2 - 2*vs^2)
        };

        PETSC_DEVICE static inline Unpacked
        unpack(const PylithInt aOff[],
               const PylithInt aOff_x[],
               const PylithScalar a[],
               const PylithScalar a_t[],
               const PylithScalar a_x[]) {
            const PylithScalar density = a[aOff[DENSITY]];
            const PylithScalar vs      = a[aOff[VS]];
            const PylithScalar vp      = a[aOff[VP]];
            return {
                density,
                density * vs * vs,                    // mu
                density * (vp*vp - 2.0*vs*vs)         // lambda
            };
        }
    };

    // ----------------------------------------------------------
    // Kernels
    // ----------------------------------------------------------

    // f0: body force (buoyancy, gravity)
    // f0_i = -rho * g_i
    template <PylithInt dim>
    PETSC_DEVICE static inline void
    f0(const PylithInt dim_,
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
       PylithScalar f0[]) {

        // constants[0..dim-1] = gravitational acceleration vector
        const auto aux = Auxiliary::unpack(aOff, aOff_x, a, a_t, a_x);

        for (PylithInt i = 0; i < dim; i++)
            f0[i] = -aux.density * constants[i];
    }

    // f1: stress tensor (momentum balance flux term)
    // f1_ij = sigma_ij = lambda * tr(eps) * delta_ij + 2*mu * eps_ij
    template <PylithInt dim>
    PETSC_DEVICE static inline void
    f1(const PylithInt dim_,
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
       PylithScalar f1[]) {

        const auto sol = Solution::unpack<dim>(sOff, sOff_x, s, s_t, s_x);
        const auto aux = Auxiliary::unpack(aOff, aOff_x, a, a_t, a_x);

        const PylithScalar tr = sol.displacement.trace();

        for (PylithInt i = 0; i < dim; i++)
        for (PylithInt j = 0; j < dim; j++)
            f1[i*dim+j] = aux.lambda * tr * (i==j)
                         + 2.0 * aux.mu * sol.displacement.symGrad(i,j);
    }

    // Jf3_uu: Jacobian of f1 w.r.t. displacement gradient
    // Jf3[i,j,k,l] = lambda*delta_ij*delta_kl + mu*(delta_ik*delta_jl + delta_il*delta_jk)
    template <PylithInt dim>
    PETSC_DEVICE static inline void
    Jf3_uu(const PylithInt dim_,
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
            const PylithReal s_tShift,
            const PylithReal x[],
            const PylithReal n[],
            const PylithInt numConstants,
            const PylithScalar constants[],
            PylithScalar Jf3[]) {

        const auto aux = Auxiliary::unpack(aOff, aOff_x, a, a_t, a_x);

        // Written directly — no C[dim][dim][dim][dim] ever materialized
        for (PylithInt i = 0; i < dim; i++)
        for (PylithInt j = 0; j < dim; j++)
        for (PylithInt k = 0; k < dim; k++)
        for (PylithInt l = 0; l < dim; l++)
            Jf3[((i*dim+j)*dim+k)*dim+l] =
                aux.lambda * (i==j) * (k==l)
                + aux.mu * ((i==k)*(j==l) + (i==l)*(j==k));
    }

    // ----------------------------------------------------------
    // Dynamic (inertial) variants
    // ----------------------------------------------------------

    // f0_inertial: rho * a_i  (inertial term, added to body force)
    // f0_i = rho * s_tt_i - rho * g_i
    template <PylithInt dim>
    PETSC_DEVICE static inline void
    f0_inertial(const PylithInt dim_,
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
                PylithScalar f0[]) {

        const auto sol = Solution::unpack<dim>(sOff, sOff_x, s, s_t, s_x);
        const auto aux = Auxiliary::unpack(aOff, aOff_x, a, a_t, a_x);

        for (PylithInt i = 0; i < dim; i++)
            f0[i] = aux.density * sol.displacement.dot(i)
                   - aux.density * constants[i];
    }

    // Jf0_uu_inertial: d(f0_inertial)/d(s_t) * s_tShift
    // Jf0_ij = rho * s_tShift * delta_ij
    template <PylithInt dim>
    PETSC_DEVICE static inline void
    Jf0_uu_inertial(const PylithInt dim_,
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
                    const PylithReal s_tShift,
                    const PylithReal x[],
                    const PylithReal n[],
                    const PylithInt numConstants,
                    const PylithScalar constants[],
                    PylithScalar Jf0[]) {

        const auto aux = Auxiliary::unpack(aOff, aOff_x, a, a_t, a_x);

        for (PylithInt i = 0; i < dim; i++)
        for (PylithInt j = 0; j < dim; j++)
            Jf0[i*dim+j] = aux.density * s_tShift * (i==j);
    }

    // ----------------------------------------------------------
    // Boundary kernels
    // ----------------------------------------------------------

    // Traction boundary condition:
    // f0_i = n_j * sigma_ij  (natural BC — outward traction)
    template <PylithInt dim>
    PETSC_DEVICE static inline void
    f0_traction(const PylithInt dim_,
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
                PylithScalar f0[]) {

        // Traction stored in aux as a dim-vector
        // aOff[0] = traction vector field
        for (PylithInt i = 0; i < dim; i++)
            f0[i] -= a[aOff[0] + i];   // sign: traction opposes residual
    }

    // ----------------------------------------------------------
    // Registration helper
    // ----------------------------------------------------------

    // Returns the correct function pointers for a given dim at runtime,
    // for registration with PetscDS / PyLith integrators.
    // Call once during problem setup — not in the hot loop.

    template <PylithInt dim>
    struct Kernels {
        static constexpr auto f0_body   = &IsotropicLinearElasticity::f0<dim>;
        static constexpr auto f1_stress = &IsotropicLinearElasticity::f1<dim>;
        static constexpr auto Jf3       = &IsotropicLinearElasticity::Jf3_uu<dim>;
        static constexpr auto f0_inert  = &IsotropicLinearElasticity::f0_inertial<dim>;
        static constexpr auto Jf0_inert = &IsotropicLinearElasticity::Jf0_uu_inertial<dim>;
    };
};


// Setup: branch on dim once, never again
if (spaceDim == 2) {
    using K = IsotropicLinearElasticity::Kernels<2>;
    PetscDSSetResidual(ds, 0, K::f0_body,   K::f1_stress);
    PetscDSSetJacobian(ds, 0, 0, nullptr, nullptr, nullptr, K::Jf3);
} else {
    using K = IsotropicLinearElasticity::Kernels<3>;
    PetscDSSetResidual(ds, 0, K::f0_body,   K::f1_stress);
    PetscDSSetJacobian(ds, 0, 0, nullptr, nullptr, nullptr, K::Jf3);
}






// Presence flags — combine with | to form a field mask
enum AuxFlags : PylithInt {
    AUX_NONE        = 0,
    AUX_BODY_FORCE  = 1 << 0,
    AUX_GRAVITY     = 1 << 1,
    AUX_REFERENCE   = 1 << 2,   // reference stress/strain
    // extensible...
};

template <AuxFlags flags>
struct AuxiliaryLayout {

    // Is a given flag present?
    static constexpr bool has(AuxFlags f) { return (flags & f) != 0; }

    // Compile-time field index: count how many required fields
    // precede this optional one, plus how many optional fields
    // before it are enabled. Required fields always come first.
    enum Fields : PylithInt {
        // Required fields — always present, always at fixed indices
        DENSITY  = 0,
        VS       = 1,
        VP       = 2,

        // Optional fields — index depends on which precede them
        BODY_FORCE = 3,
        GRAVITY    = 3 + (has(AUX_BODY_FORCE) ? 1 : 0),
        REFERENCE  = 3 + (has(AUX_BODY_FORCE) ? 1 : 0)
                       + (has(AUX_GRAVITY)     ? 1 : 0),

        NUM_FIELDS = 3 + (has(AUX_BODY_FORCE) ? 1 : 0)
                       + (has(AUX_GRAVITY)     ? 1 : 0)
                       + (has(AUX_REFERENCE)   ? 1 : 0),
    };

    // ----------------------------------------------------------
    // Unpacked: only contains members for present fields
    // Absent fields take zero space — no dead pointers
    // ----------------------------------------------------------

    // Helper: conditionally include a member
    // Empty base class optimization ensures zero size when unused
    template <bool present, typename T>
    struct OptionalMember {};

    template <typename T>
    struct OptionalMember<true, T> { T value; };

    template <PylithInt dim>
    struct Unpacked {
        // Required — always present
        PylithScalar density;
        PylithScalar mu;
        PylithScalar lambda;

        // Optional — zero size if absent
        [[no_unique_address]]
        OptionalMember<has(AUX_BODY_FORCE), VectorField<dim>> body_force;

        [[no_unique_address]]
        OptionalMember<has(AUX_GRAVITY),    VectorField<dim>> gravity;

        [[no_unique_address]]
        OptionalMember<has(AUX_REFERENCE),  TensorField<dim>> reference_stress;

        // Type-safe accessor — compile error if field not present
        template <AuxFlags F>
        PETSC_DEVICE inline auto& get() {
            static_assert(has(F), "Auxiliary field not present in this layout");
            if constexpr (F == AUX_BODY_FORCE) return body_force.value;
            if constexpr (F == AUX_GRAVITY)    return gravity.value;
            if constexpr (F == AUX_REFERENCE)  return reference_stress.value;
        }
    };

    // ----------------------------------------------------------
    // unpack(): builds Unpacked from raw arrays
    // Optional field unpacking is compiled away when absent
    // ----------------------------------------------------------

    template <PylithInt dim>
    PETSC_DEVICE static inline Unpacked<dim>
    unpack(const PylithInt  aOff[],
           const PylithInt  aOff_x[],
           const PylithScalar a[],
           const PylithScalar a_t[],
           const PylithScalar a_x[]) {

        const PylithScalar density = a[aOff[DENSITY]];
        const PylithScalar vs      = a[aOff[VS]];
        const PylithScalar vp      = a[aOff[VP]];

        Unpacked<dim> u;
        u.density = density;
        u.mu      = density * vs * vs;
        u.lambda  = density * (vp*vp - 2.0*vs*vs);

        if constexpr (has(AUX_BODY_FORCE))
            u.body_force.value = { &a[aOff[BODY_FORCE]], nullptr, nullptr };

        if constexpr (has(AUX_GRAVITY))
            u.gravity.value    = { &a[aOff[GRAVITY]], nullptr, nullptr };

        if constexpr (has(AUX_REFERENCE))
            u.reference_stress.value = { &a[aOff[REFERENCE]] };

        return u;
    }
};


template <PylithInt dim, AuxFlags flags>
struct ElasticityKernels {
    using AL = AuxiliaryLayout<flags>;
    using SL = IsotropicLinearElasticity::Solution;

    // f0: body terms — gravity and/or body force, whichever are present
    PETSC_DEVICE static inline void
    f0(const PylithInt dim_,
       const PylithInt numS,       const PylithInt numA,
       const PylithInt sOff[],     const PylithInt sOff_x[],
       const PylithScalar s[],     const PylithScalar s_t[],
       const PylithScalar s_x[],
       const PylithInt aOff[],     const PylithInt aOff_x[],
       const PylithScalar a[],     const PylithScalar a_t[],
       const PylithScalar a_x[],
       const PylithReal t,         const PylithReal x[],
       const PylithReal n[],
       const PylithInt numConstants, const PylithScalar constants[],
       PylithScalar f0[]) {

        const auto aux = AL::template unpack<dim>(aOff, aOff_x, a, a_t, a_x);

        // Initialize to zero
        for (PylithInt i = 0; i < dim; i++) f0[i] = 0.0;

        // Gravity: f0_i += rho * g_i
        // Compiled away entirely if AUX_GRAVITY not set
        if constexpr (AL::has(AUX_GRAVITY)) {
            const auto& g = aux.template get<AUX_GRAVITY>();
            for (PylithInt i = 0; i < dim; i++)
                f0[i] -= aux.density * g(i);
        }

        // Body force: f0_i += b_i
        // Compiled away entirely if AUX_BODY_FORCE not set
        if constexpr (AL::has(AUX_BODY_FORCE)) {
            const auto& b = aux.template get<AUX_BODY_FORCE>();
            for (PylithInt i = 0; i < dim; i++)
                f0[i] -= b(i);
        }
    }

    // f1: stress — optionally adds reference stress
    PETSC_DEVICE static inline void
    f1(const PylithInt dim_,
       const PylithInt numS,       const PylithInt numA,
       const PylithInt sOff[],     const PylithInt sOff_x[],
       const PylithScalar s[],     const PylithScalar s_t[],
       const PylithScalar s_x[],
       const PylithInt aOff[],     const PylithInt aOff_x[],
       const PylithScalar a[],     const PylithScalar a_t[],
       const PylithScalar a_x[],
       const PylithReal t,         const PylithReal x[],
       const PylithReal n[],
       const PylithInt numConstants, const PylithScalar constants[],
       PylithScalar f1[]) {

        const auto sol = SL::template unpack<dim>(sOff, sOff_x, s, s_t, s_x);
        const auto aux = AL::template unpack<dim>(aOff, aOff_x, a, a_t, a_x);

        const PylithScalar tr = sol.displacement.trace();

        for (PylithInt i = 0; i < dim; i++)
        for (PylithInt j = 0; j < dim; j++)
            f1[i*dim+j] = aux.lambda * tr * (i==j)
                         + 2.0 * aux.mu * sol.displacement.symGrad(i,j);

        // Reference stress subtracted from residual
        if constexpr (AL::has(AUX_REFERENCE)) {
            const auto& s0 = aux.template get<AUX_REFERENCE>();
            for (PylithInt i = 0; i < dim; i++)
            for (PylithInt j = 0; j < dim; j++)
                f1[i*dim+j] -= s0(i,j);
        }
    }
};



// Forward declaration of the registration worker
template <PylithInt dim, AuxFlags flags>
void registerKernels(PetscDS ds, PylithInt field, bool dynamic);

// Runtime dispatch: branch on dim and flags once at setup
void setupElasticity(PetscDS ds, PylithInt dim, PylithInt field,
                     bool hasGravity, bool hasBodyForce,
                     bool hasReference, bool dynamic) {

    const AuxFlags flags = static_cast<AuxFlags>(
          (hasGravity    ? AUX_GRAVITY     : AUX_NONE)
        | (hasBodyForce  ? AUX_BODY_FORCE  : AUX_NONE)
        | (hasReference  ? AUX_REFERENCE   : AUX_NONE));

    // Expand all combinations — compiler emits only reachable paths
    #define REGISTER(D, F) \
        if (dim == D && flags == F) { registerKernels<D,F>(ds, field, dynamic); return; }

    REGISTER(2, AUX_NONE)
    REGISTER(2, AUX_GRAVITY)
    REGISTER(2, AUX_BODY_FORCE)
    REGISTER(2, AUX_GRAVITY | AUX_BODY_FORCE)
    REGISTER(2, AUX_REFERENCE)
    REGISTER(2, AUX_GRAVITY | AUX_REFERENCE)
    REGISTER(2, AUX_BODY_FORCE | AUX_REFERENCE)
    REGISTER(2, AUX_GRAVITY | AUX_BODY_FORCE | AUX_REFERENCE)
    REGISTER(3, AUX_NONE)
    REGISTER(3, AUX_GRAVITY)
    REGISTER(3, AUX_BODY_FORCE)
    REGISTER(3, AUX_GRAVITY | AUX_BODY_FORCE)
    REGISTER(3, AUX_REFERENCE)
    REGISTER(3, AUX_GRAVITY | AUX_REFERENCE)
    REGISTER(3, AUX_BODY_FORCE | AUX_REFERENCE)
    REGISTER(3, AUX_GRAVITY | AUX_BODY_FORCE | AUX_REFERENCE)
    #undef REGISTER
}

template <PylithInt dim, AuxFlags flags>
void registerKernels(PetscDS ds, PylithInt field, bool dynamic) {
    using K = ElasticityKernels<dim, flags>;
    PetscDSSetResidual(ds, field, K::f0, K::f1);
    PetscDSSetJacobian(ds, field, field,
                       nullptr, nullptr, nullptr,
                       K::Jf3_uu);
    if (dynamic)
        PetscDSSetJacobian(ds, field, field,
                           K::Jf0_uu_inertial,
                           nullptr, nullptr, nullptr);
}
