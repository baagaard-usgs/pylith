# GreensFns

% WARNING: Do not edit; this is a generated file!
Full name: `pylith.problems.GreensFns`

Static Green's function problem type with each Green's function corresponding to a fault slip impulses.

Implements `Problem`.

## Pyre Facilities

* `bc`: Boundary conditions.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, bc
* `defaults`: Default options for problem.
  - **current value**: 'problem_defaults', from {default}
  - **configurable as**: problem_defaults, defaults
* `gravity_field`: Database used for gravity field.
  - **current value**: 'nullcomponent', from {default}
  - **configurable as**: nullcomponent, gravity_field
* `interfaces`: Interior surfaces with constraints or constitutive models.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, interfaces
* `materials`: Materials in problem.
  - **current value**: 'homogeneous', from {default}
  - **configurable as**: homogeneous, materials
* `normalizer`: Nondimensionalizer for problem.
  - **current value**: 'nondimelasticquasistatic', from {default}
  - **configurable as**: nondimelasticquasistatic, normalizer
* `progress_monitor`: Simple progress monitor via text file.
  - **current value**: 'progressmonitorstep', from {default}
  - **configurable as**: progressmonitorstep, progress_monitor
* `solution`: Solution field for problem.
  - **current value**: 'solution', from {default}
  - **configurable as**: solution
* `solution_observers`: Observers (e.g., output) for solution.
  - **current value**: 'singlesolnobserver', from {default}
  - **configurable as**: singlesolnobserver, solution_observers

## Pyre Properties

* `formulation`=\<str\>: Formulation for equations.
  - **default value**: 'quasistatic'
  - **current value**: 'quasistatic', from {default}
  - **validator**: (in ['quasistatic', 'dynamic', 'dynamic_imex'])
* `label`=\<str\>: Name of label identifier for fault surface on which to impose impulses.
  - **default value**: 'fault'
  - **current value**: 'fault', from {default}
* `label_value`=\<int\>: Value of label identifier for fault surface on which to impose impulses.
  - **default value**: 1
  - **current value**: 1, from {default}
* `solver`=\<str\>: Type of solver to use ['linear', 'nonlinear'].
  - **default value**: 'linear'
  - **current value**: 'linear', from {default}
  - **validator**: (in ['linear', 'nonlinear'])

## Example

Example of setting `GreensFns` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp]
problem = pylith.problems.GreensFns

[pylithapp.greensfns]
label = fault
label_value = 1

interfaces = [fault]
interfaces.fault = pylith.faults.FaultCohesiveImpulses

[pylithapp.greensfns.interfaces.fault]
label = fault
label_value = 20

# Impulses for left-lateral slip (dof=1)
impulse_dof = [1]
:::
