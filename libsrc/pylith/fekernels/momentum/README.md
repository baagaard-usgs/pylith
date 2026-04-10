# Momentum equation

## Terms

- body force
- divergence of stress
- inertia

## SolutionLayout

- Required subfields
  - `displacement` VectorField

## AuxiliaryLayout

- Required subfields
  - `density` ScalarField
- Optional subfields
  - `body_force` VectorField
  - `gravitational_acceleration` VectorField

## StrainModel

- Required methods
  - `compute(Matrix<dim>& strain, const SolutionLayout& solution)`

## StressModel

- Required methods
  - `cauchyStress(Matrix<dim>& stress, const Matrix<dim>& strain, const AuxiliaryLayout& auxiliary)`
  - `cauchyStressTangent(pylith::scalar Jf3[], const AuxiliaryLayout& auxiliary, const pylith::scalar sign)`
