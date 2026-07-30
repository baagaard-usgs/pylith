(sec-user-governing-eqns-formulation)=
# A Common Finite-Element Formulation

Every governing equation in PyLith — elasticity, incompressible elasticity, and poroelasticity — is solved using the same finite-element formulation.
We cast each governing equation into a common weak form, $F(t,s,\dot{s}) = G(t,s)$, and solve it using the PETSc finite-element and time-stepping infrastructure.
What changes from one governing equation to another is *not* the machinery, but two things: the set of solution subfields (for example, displacement alone for elasticity, versus displacement, pressure, and volumetric strain for poroelasticity) and the pointwise functions that define the terms in the weak form.
Because the formulation is shared, the residual and Jacobian for each governing equation in the sections that follow are expressed using the same notation introduced here.

Within the PETSc solver framework, we want to solve a system of partial differential equations in which the weak form can be expressed as $F(t,s,\dot{s}) = G(t,s)$, $s(t_0) = s_0$, where $F$ and $G$ are vector functions, $t$ is time, and $s$ is the solution vector.

% ... existing content unchanged: the general weak form (Eq. for the integrals),
% the identification of G as the RHS function and F as the LHS (or I) function,
% the pointwise-function decomposition into f0, f1, g0, g1, the Jacobian
% structure J0-J3, and the PETSc TS notes on explicit, implicit, and
% implicit-explicit time stepping ...

:::{note}
The pointwise functions $\vec{f}_0$, $\boldsymbol{f}_1$, $\vec{g}_0$, $\boldsymbol{g}_1$ and the Jacobian pointwise functions $J_0$ through $J_3$ introduced above are the vocabulary used throughout the rest of this chapter.
Each governing equation is presented in the same order: we state the physics, give the strong form, derive the weak form, and then identify the pointwise functions using this notation.
This makes it straightforward to see how each governing equation maps onto the shared formulation and to compare terms across equations.
:::

