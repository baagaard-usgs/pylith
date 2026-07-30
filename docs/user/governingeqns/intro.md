(sec-user-governing-eqns)=
# Governing Equations

This chapter describes the governing equations that PyLith solves and how we cast each one into a common finite-element formulation.
You can read this chapter at two levels of detail.
If you want to understand the physics and know which problems PyLith can solve, read the {ref}`overview <sec-user-governing-eqns-overview>` and the introductory discussion at the beginning of each governing equation.
If you want to verify exactly what the code solves or extend it, continue into the weak forms and pointwise functions that follow.

(sec-user-governing-eqns-overview)=
## Overview and Scope

PyLith simulates crustal deformation across a wide range of spatial and temporal scales.
The governing equations describe the balance of momentum in a deforming solid and, for poroelasticity, the coupled flow of fluid through that solid.
In all cases we solve the equations using the finite-element method, casting each governing equation into the common form described in {ref}`sec-user-governing-eqns-formulation`.

### Implemented governing equations

PyLith currently implements the three governing equations listed in {numref}`tab-governing-eqns-implemented`.
Each governing equation supports one or more *formulations* depending on whether inertia is included (dynamic), neglected (quasistatic), or the problem is time independent (static).

```{table} Governing equations implemented in PyLith.
:name: tab-governing-eqns-implemented

| Governing equation | Formulations | Typical use |
| :--- | :--- | :--- |
| Elasticity | static, quasistatic, dynamic | Interseismic and coseismic deformation; seismic wave propagation |
| Incompressible elasticity | quasistatic | Estimating a lithostatic (gravitational) stress state with minimal volumetric deformation |
| Poroelasticity | quasistatic, dynamic | Coupled deformation and fluid flow (for example, magma reservoirs, fluid-driven faulting) |
```

For elasticity we support a variety of bulk rheologies, including linear elastic, linear Maxwell viscoelastic, generalized Maxwell viscoelastic, and power-law viscoelastic behavior.
Poroelasticity uses an isotropic, linear poroelastic rheology and includes an optional porosity state variable that evolves in time.

### Faults as interior interfaces

PyLith treats faults as interior interfaces embedded in the domain.
We can prescribe slip on these interfaces for the elasticity and poroelasticity equations.
Prescribed (kinematic) slip is currently supported; spontaneous rupture (fault friction) is planned for a future release (see {ref}`sec-user-governing-eqns` limitations and the Development Plan).
The mathematical treatment of fault interfaces is included with each governing equation that supports them, and the physics user interface is described in {ref}`sec-user-physics-faults`.

### Mathematical notation

We use the notation in {numref}`tab-governing-eqns-notation` throughout this chapter.

```{table} Mathematical notation for governing equations.
:name: tab-governing-eqns-notation

| Symbol | Description |
| :--- | :--- |
| $\vec{a}$ | Vector field $a$ |
| $\mathbf{a}$ | Second-order tensor field $a$ |
| $\vec{u}$ | Displacement vector field |
| $\vec{v}$ | Velocity vector field |
| $\vec{d}$ | Fault slip vector field |
| $\vec{f}$ | Body force vector field |
| $\vec{\tau}$ | Traction vector field |
| $\boldsymbol{\sigma}$ | Stress tensor field |
| $\vec{n}$ | Normal vector field |
| $\rho$ | Mass density scalar field |
```

:::{seealso}
The nondimensionalization of each governing equation is discussed within the corresponding section.
The user interface for selecting governing equations, rheologies, and their parameters is described in {ref}`sec-user-physics-materials`.
:::

