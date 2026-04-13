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



REGISTRY

```
ProblemConfig (runtime)
    │
    ├─ solution_layout() ──► SolutionLayout enum
    ├─ aux_layout()      ──► AuxLayout enum
    │
    └─► KernelKey { dim, solution, aux, strain, stress }
              │
              │  hash lookup (O(1))
              ▼
        KernelRegistry::map_
              │
              │  populated at startup by walking ValidConfigs
              │  register_one<Config<Dim,S,T>>() for each entry
              ▼
        &kernel<Dim, StrainModel, StressModel>   (function pointer)
              │
              ▼
        fn(data, n)
```


```
#include <unordered_map>
#include <stdexcept>
#include <tuple>
#include <cstdint>
#include <string>

// ============================================================
// 1. DIMENSIONS
// ============================================================

struct Dim2 { static constexpr int value = 2; };
struct Dim3 { static constexpr int value = 3; };

// ============================================================
// 2. LAYOUT ENUMS
// ============================================================

enum class SolutionLayout : uint32_t {
    DEFAULT           = 0,
    FAULT             = 1 << 0,
    INERTIA           = 1 << 1,
    FAULT_AND_INERTIA = FAULT | INERTIA,
};

enum class AuxLayout : uint32_t {
    DEFAULT                = 0,
    BODY_FORCE             = 1 << 0,
    GRAVITY                = 1 << 1,
    BODY_FORCE_AND_GRAVITY = BODY_FORCE | GRAVITY,
};

// ============================================================
// 3. STRAIN AND STRESS MODEL TEMPLATES
// ============================================================

template<typename Dim, SolutionLayout S> struct Infinitesimal {};
template<typename Dim, SolutionLayout S> struct Finite {};

template<typename Dim, AuxLayout A> struct IsotropicLinear {};
template<typename Dim, AuxLayout A> struct IsotropicMaxwell {};
template<typename Dim, AuxLayout A> struct IsotropicPowerLaw {};

// ============================================================
// 4. RUNTIME ENUMS AND PROBLEM CONFIG
// ============================================================

enum class StrainType  { Infinitesimal, Finite };
enum class StressType  { IsotropicLinear, IsotropicMaxwell, IsotropicPowerLaw };

struct ProblemConfig {
    int        dim;
    bool       has_fault;
    bool       has_inertia;
    bool       has_body_force;
    bool       has_gravity;
    StrainType strain;
    StressType stress;
};

SolutionLayout solution_layout(const ProblemConfig& cfg) {
    uint32_t v = 0;
    if (cfg.has_fault)   v |= uint32_t(SolutionLayout::FAULT);
    if (cfg.has_inertia) v |= uint32_t(SolutionLayout::INERTIA);
    return SolutionLayout(v);
}

AuxLayout aux_layout(const ProblemConfig& cfg) {
    uint32_t v = 0;
    if (cfg.has_body_force) v |= uint32_t(AuxLayout::BODY_FORCE);
    if (cfg.has_gravity)    v |= uint32_t(AuxLayout::GRAVITY);
    return AuxLayout(v);
}

// ============================================================
// 5. RUNTIME KEY
// ============================================================

struct KernelKey {
    int            dim;
    SolutionLayout solution;
    AuxLayout      aux;
    StrainType     strain;
    StressType     stress;

    bool operator==(const KernelKey&) const = default;
};

struct KernelKeyHash {
    size_t operator()(const KernelKey& k) const {
        auto h = [](size_t seed, auto v) {
            return seed ^ (std::hash<size_t>{}(size_t(v)) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
        };
        size_t s = std::hash<int>{}(k.dim);
        s = h(s, k.solution);
        s = h(s, k.aux);
        s = h(s, k.strain);
        s = h(s, k.stress);
        return s;
    }
};

// ============================================================
// 6. MODEL TRAITS — extract runtime identity from template types
// ============================================================

template<typename T> struct model_traits;

template<typename Dim, SolutionLayout S>
struct model_traits<Infinitesimal<Dim, S>> {
    static constexpr StrainType     strain   = StrainType::Infinitesimal;
    static constexpr SolutionLayout solution = S;
};

template<typename Dim, SolutionLayout S>
struct model_traits<Finite<Dim, S>> {
    static constexpr StrainType     strain   = StrainType::Finite;
    static constexpr SolutionLayout solution = S;
};

template<typename Dim, AuxLayout A>
struct model_traits<IsotropicLinear<Dim, A>> {
    static constexpr StressType stress = StressType::IsotropicLinear;
    static constexpr AuxLayout  aux    = A;
};

template<typename Dim, AuxLayout A>
struct model_traits<IsotropicMaxwell<Dim, A>> {
    static constexpr StressType stress = StressType::IsotropicMaxwell;
    static constexpr AuxLayout  aux    = A;
};

template<typename Dim, AuxLayout A>
struct model_traits<IsotropicPowerLaw<Dim, A>> {
    static constexpr StressType stress = StressType::IsotropicPowerLaw;
    static constexpr AuxLayout  aux    = A;
};

// ============================================================
// 7. CONFIG TYPE AND KEY EXTRACTION
// ============================================================

template<typename Dim, typename StrainModel, typename StressModel>
struct Config {};

template<typename Dim, typename StrainModel, typename StressModel>
constexpr KernelKey key_for() {
    return {
        Dim::value,
        model_traits<StrainModel>::solution,
        model_traits<StressModel>::aux,
        model_traits<StrainModel>::strain,
        model_traits<StressModel>::stress,
    };
}

// ============================================================
// 8. KERNEL
// ============================================================

using KernelFn = void(*)(float*, size_t);

template<typename Dim, typename StrainModel, typename StressModel>
void kernel(float* data, size_t n) {
    // Compile-time branches — dead code is eliminated
    if constexpr (std::is_same_v<Dim, Dim2>)                              { /* 2D setup   */ }
    if constexpr (std::is_same_v<Dim, Dim3>)                              { /* 3D setup   */ }
    if constexpr (model_traits<StrainModel>::strain == StrainType::Finite) { /* finite     */ }
    if constexpr (model_traits<StressModel>::stress == StressType::IsotropicMaxwell) { /* maxwell */ }
    if constexpr ((uint32_t(model_traits<StressModel>::aux) &
                   uint32_t(AuxLayout::GRAVITY)) != 0)                    { /* gravity    */ }
}

// ============================================================
// 9. VALID CONFIGS — the physics contract
// ============================================================

using ValidConfigs = std::tuple
    // 2D — Infinitesimal strain
    Config<Dim2, Infinitesimal<Dim2, SolutionLayout::DEFAULT>,           IsotropicLinear<Dim2,   AuxLayout::DEFAULT>>,
    Config<Dim2, Infinitesimal<Dim2, SolutionLayout::DEFAULT>,           IsotropicLinear<Dim2,   AuxLayout::BODY_FORCE>>,
    Config<Dim2, Infinitesimal<Dim2, SolutionLayout::DEFAULT>,           IsotropicLinear<Dim2,   AuxLayout::GRAVITY>>,
    Config<Dim2, Infinitesimal<Dim2, SolutionLayout::DEFAULT>,           IsotropicLinear<Dim2,   AuxLayout::BODY_FORCE_AND_GRAVITY>>,
    Config<Dim2, Infinitesimal<Dim2, SolutionLayout::FAULT>,             IsotropicLinear<Dim2,   AuxLayout::DEFAULT>>,
    Config<Dim2, Infinitesimal<Dim2, SolutionLayout::FAULT>,             IsotropicLinear<Dim2,   AuxLayout::GRAVITY>>,
    Config<Dim2, Infinitesimal<Dim2, SolutionLayout::DEFAULT>,           IsotropicMaxwell<Dim2,  AuxLayout::DEFAULT>>,
    Config<Dim2, Infinitesimal<Dim2, SolutionLayout::DEFAULT>,           IsotropicMaxwell<Dim2,  AuxLayout::GRAVITY>>,
    Config<Dim2, Infinitesimal<Dim2, SolutionLayout::DEFAULT>,           IsotropicPowerLaw<Dim2, AuxLayout::DEFAULT>>,
    Config<Dim2, Finite<Dim2,        SolutionLayout::DEFAULT>,           IsotropicLinear<Dim2,   AuxLayout::DEFAULT>>,
    Config<Dim2, Finite<Dim2,        SolutionLayout::DEFAULT>,           IsotropicLinear<Dim2,   AuxLayout::GRAVITY>>,

    // 3D — Infinitesimal strain
    Config<Dim3, Infinitesimal<Dim3, SolutionLayout::DEFAULT>,           IsotropicLinear<Dim3,   AuxLayout::DEFAULT>>,
    Config<Dim3, Infinitesimal<Dim3, SolutionLayout::DEFAULT>,           IsotropicLinear<Dim3,   AuxLayout::BODY_FORCE>>,
    Config<Dim3, Infinitesimal<Dim3, SolutionLayout::DEFAULT>,           IsotropicLinear<Dim3,   AuxLayout::GRAVITY>>,
    Config<Dim3, Infinitesimal<Dim3, SolutionLayout::DEFAULT>,           IsotropicLinear<Dim3,   AuxLayout::BODY_FORCE_AND_GRAVITY>>,
    Config<Dim3, Infinitesimal<Dim3, SolutionLayout::FAULT>,             IsotropicLinear<Dim3,   AuxLayout::DEFAULT>>,
    Config<Dim3, Infinitesimal<Dim3, SolutionLayout::FAULT>,             IsotropicLinear<Dim3,   AuxLayout::GRAVITY>>,
    Config<Dim3, Infinitesimal<Dim3, SolutionLayout::DEFAULT>,           IsotropicMaxwell<Dim3,  AuxLayout::DEFAULT>>,
    Config<Dim3, Infinitesimal<Dim3, SolutionLayout::DEFAULT>,           IsotropicMaxwell<Dim3,  AuxLayout::GRAVITY>>,
    Config<Dim3, Infinitesimal<Dim3, SolutionLayout::DEFAULT>,           IsotropicPowerLaw<Dim3, AuxLayout::DEFAULT>>,
    Config<Dim3, Finite<Dim3,        SolutionLayout::DEFAULT>,           IsotropicLinear<Dim3,   AuxLayout::DEFAULT>>,
    Config<Dim3, Finite<Dim3,        SolutionLayout::DEFAULT>,           IsotropicLinear<Dim3,   AuxLayout::GRAVITY>>
>;

// ============================================================
// 10. KERNEL REGISTRY
// ============================================================

class KernelRegistry {
    std::unordered_map<KernelKey, KernelFn, KernelKeyHash> map_;

public:
    KernelRegistry() {
        register_all(std::make_index_sequence<std::tuple_size_v<ValidConfigs>>{});
    }

    KernelFn get(const ProblemConfig& cfg) const {
        KernelKey key {
            cfg.dim,
            solution_layout(cfg),
            aux_layout(cfg),
            cfg.strain,
            cfg.stress,
        };
        auto it = map_.find(key);
        if (it == map_.end())
            throw std::invalid_argument(unsupported_msg(key));
        return it->second;
    }

private:
    template<size_t... Is>
    void register_all(std::index_sequence<Is...>) {
        (register_one<std::tuple_element_t<Is, ValidConfigs>>(), ...);
    }

    template<typename Cfg>
    void register_one();                           // defined below via partial specialization
};

template<typename Dim, typename StrainModel, typename StressModel>
void register_config(std::unordered_map<KernelKey, KernelFn, KernelKeyHash>& map) {
    map.emplace(key_for<Dim, StrainModel, StressModel>(),
                &kernel<Dim, StrainModel, StressModel>);
}

template<> template<typename Dim, typename S, typename T>
void KernelRegistry::register_one<Config<Dim, S, T>>() {
    register_config<Dim, S, T>(map_);
}

std::string unsupported_msg(const KernelKey& k) {
    return "Unsupported physics combination: dim=" + std::to_string(k.dim)
         + " solution=" + std::to_string(uint32_t(k.solution))
         + " aux="      + std::to_string(uint32_t(k.aux))
         + " strain="   + std::to_string(int(k.strain))
         + " stress="   + std::to_string(int(k.stress));
}

// ============================================================
// 11. USAGE
// ============================================================

static const KernelRegistry registry;

void solve(const ProblemConfig& cfg, float* data, size_t n) {
    KernelFn fn = registry.get(cfg);
    fn(data, n);
}

int main() {
    ProblemConfig cfg {
        .dim           = 2,
        .has_fault     = false,
        .has_inertia   = false,
        .has_body_force = false,
        .has_gravity   = true,
        .strain        = StrainType::Infinitesimal,
        .stress        = StressType::IsotropicMaxwell,
    };

    float data[1024] = {};
    solve(cfg, data, 1024);   // dispatches to kernel<Dim2, Infinitesimal<Dim2, DEFAULT>, IsotropicMaxwell<Dim2, GRAVITY>>

    return 0;
}
```