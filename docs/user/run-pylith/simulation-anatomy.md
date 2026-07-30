(sec-user-anatomy)=
# Anatomy of a Simulation

Before you write your own parameter files, it helps to have a mental picture of what a PyLith simulation is made of.
Every simulation, whether it is a static box being stretched or a 3D subduction zone rupturing over thousands of years, is assembled from the same handful of ingredients.
This page walks through those ingredients, how they fit together, and where each one lives in the parameter files.

## The big picture

A PyLith simulation takes three things you provide and produces output files you analyze:

```{mermaid}
:config: {"theme":"neutral"}
flowchart LR
    subgraph inputs["You provide"]
        mesh["Finite-element mesh<br/>(geometry + topology)"]
        params["Parameters<br/>(.cfg files)"]
        dbs["Spatial databases<br/>(property values)"]
    end
    pylith(["PyLith"])
    subgraph outputs["PyLith produces"]
        soln["Solution fields<br/>(.h5 / .vtu)"]
        derived["Derived fields<br/>(stress, strain)"]
    end
    mesh --> pylith
    params --> pylith
    dbs --> pylith
    pylith --> soln
    pylith --> derived
```

The mesh defines *where* the problem lives.
The parameters define *what physics* to solve and *how*.
The spatial databases define *what the material properties and boundary values are* at each location.
PyLith combines these, solves the governing equations, and writes the results.

## The five ingredients of a problem

Inside the parameters, a boundary value problem is built from five kinds of components.
Think of the first four as answering "what is being modeled" and the last as "what do I want to see."

**The mesh** is the discretized geometry of your domain, imported from Gmsh, Cubit, or an ASCII file.
Each cell carries a material identifier, and groups of cell faces or vertices are labeled so you can attach boundary conditions and faults to them.
The mesh is independent of the physics: the same mesh can be used for an elastic problem or a poroelastic one.

**Materials** occupy the volume of the domain (or the area, in 2D).
Each material specifies a governing equation, such as elasticity or poroelasticity, together with a bulk rheology, such as isotropic linear elasticity or Maxwell viscoelasticity.
A material knows which cells it applies to through its label value, which matches the material identifier in the mesh.

**Boundary conditions** act on the external surfaces of the domain.
Dirichlet conditions prescribe values (for example, displacement), Neumann conditions prescribe tractions, and absorbing conditions prevent waves from reflecting off the edges in dynamic simulations.
Each one attaches to a labeled boundary in the mesh.

**Fault interfaces** are internal surfaces where the two sides can move relative to each other.
PyLith inserts zero-thickness cohesive cells along these surfaces so it can control the jump in the solution across them.
Faults are optional; many simulations have none.

**Observers** decide what output is written and how often.
There are observers for the solution over the whole domain, over a boundary, or at discrete points, plus observers attached to each material and fault.
Observers do not affect the answer; they only record it.

```{mermaid}
:config: {"theme":"neutral"}
flowchart TB
    problem["Problem<br/>(TimeDependent or GreensFns)"]
    problem --> soln["Solution field<br/>(the unknowns you solve for)"]
    problem --> mats["Materials"]
    problem --> bcs["Boundary conditions"]
    problem --> faults["Fault interfaces<br/>(optional)"]
    problem --> obs["Solution observers"]
    mats --> mobs["Material observers"]
    faults --> fobs["Fault observers"]
```

## The solution field ties it together

The one component that connects all of the physics is the **solution field**: the set of unknowns PyLith actually solves for.
For simple elasticity the solution is just displacement.
Add a fault and you also solve for a Lagrange multiplier that represents the fault tractions.
Switch to poroelasticity and the solution grows to include fluid pressure and volumetric strain.

Choosing the right solution field is one of the most important decisions when defining a problem, because it must contain every unknown that your materials and boundary conditions require.
PyLith provides ready-made containers for common combinations (for example, `SolnDisp` for displacement alone, or `SolnDispLagrange` for displacement with a fault), so most of the time you select one rather than build it yourself.

## How the pieces map to parameter files

Each ingredient corresponds to a section in a `.cfg` file.
A skeleton for a problem with one material, two boundary conditions, and one fault looks like this:

```{code-block} cfg
[pylithapp.problem]
# The unknowns we solve for.
solution = pylith.problems.SolnDispLagrange

# Where output over the whole domain goes.
solution_observers = [domain]

[pylithapp.problem]
# One material, occupying cells with the matching label value.
materials = [elastic]

[pylithapp.problem.materials.elastic]
label_value = 1
# ... governing equation, rheology, and property database go here ...

[pylithapp.problem]
# Two external boundaries.
bc = [bc_xneg, bc_xpos]

[pylithapp.problem]
# One internal fault interface.
interfaces = [fault]
```

The property *values* (densities, moduli, slip, boundary displacements) usually do not appear in the `.cfg` file directly.
Instead, each component points to a **spatial database** that supplies values as a function of position.
This separation means you can refine your mesh or change the geometry without editing your material properties, and vice versa.

## What PyLith does with all of this

When you run `pylith`, the pieces are assembled in a fixed order.
PyLith reads the mesh, inserts cohesive cells for any faults, distributes the mesh across processes if you are running in parallel, and then hands the problem to a solver.
During the solve it queries your spatial databases to populate material properties and boundary values, integrates the governing equations over each material, applies the boundary and fault conditions, and advances the solution in time.
At each output step, the observers write the requested fields to disk.

You do not need to understand the numerical machinery to define a working simulation.
What matters at this stage is recognizing the five ingredients, knowing that the solution field must match the physics you have chosen, and understanding that properties live in spatial databases rather than in the parameter files themselves.

```{admonition} Where to go next
:class: tip
The remaining pages in this section describe each ingredient in detail: the {ref}`types of simulations <sec-user-problems>`, the {ref}`finite-element mesh <sec-user-femesh>`, and {ref}`output and observers <sec-user-output>`.
If you have not yet run a simulation, work through the box-2d examples first; seeing these ingredients in a complete, working problem makes the descriptions here much easier to follow.
```

