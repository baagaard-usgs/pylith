// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "SimpleShearPVE.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

namespace pylith {
    class _SimpleShearPVE;
} // pylith

class pylith::_SimpleShearPVE {
    static const double LENGTH_SCALE;
    static const double TIME_SCALE;
    static const double PRESSURE_SCALE;

    static const double PRESSURE0; // dimensional
    static const double XMAX; // dimensional

    // Density
    static double solid_density(const double x,
                                const double y) {
        return 2500.0;
    } // solid_density

    static double fluid_density(const double x,
                                const double y) {
        return 1000.0;
    } // fluid_density

    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    // Viscosity
    static double fluid_viscosity(const double x,
                                  const double y) {
        return 1.0e-3;
    } // vs

    static double solid_viscosity(const double x,
                                  const double y) {
        return 1.0e21;
    } // solid viscosity

    static const char* viscosity_units(void) {
        return "Pa*s";
    } // fluid_viscosity_units

    // Porosity
    static double porosity(const double x,
                           const double y) {
        return 0.02;
    } // porosity

    static const char* porosity_units(void) {
        return "none";
    } // porosity_units

    // Shear modulus
    static double shear_modulus(const double x,
                                const double y) {
        return 3.0e+10;
    } // shear_modulus

    static const char* modulus_units(void) {
        return "Pa";
    } // modulus_units

    // Drained bulk modulus
    static double drained_bulk_modulus(const double x,
                                       const double y) {
        return 8.0e+10;
    } // drained_bulk_modulus

    // Biot coefficient
    static double biot_coefficient(const double x,
                                   const double y) {
        return 1.0;
    } // biot_coefficient

    static const char* biot_coefficient_units(void) {
        return "none";
    } // biot_coefficient_units

    // Fluid modulus
    static double fluid_bulk_modulus(const double x,
                                     const double y) {
        return 1.0e+10;
    } // fluid_bulk_modulus

    // Solid modulus
    static double solid_bulk_modulus(const double x,
                                     const double y) {
        return 9.0e+10;
    } // solid_bulk_modulus

    // Permeability
    static double isotropic_permeability(const double x,
                                         const double y) {
        return 1.0e-14;
    } // isotropic_permeability

    static const char* permeability_units(void) {
        return "m*m";
    } // permeability_units

    // Strain
    static double viscous_strain_xx(const double x,
                                    const double y) {
        return 0.0;
    } // viscous_strain_xx

    static double viscous_strain_yy(const double x,
                                    const double y) {
        return 0.0;
    } // viscous_strain_yy

    static double viscous_strain_zz(const double x,
                                    const double y) {
        return 0.0;
    } // viscous_strain_zz

    static double viscous_strain_xy(const double x,
                                    const double y) {
        return -0.5;
    } // viscous_strain_xy

    static double viscous_strain_xz(const double x,
                                    const double y) {
        return 0.0;
    } // viscous_strain_xz

    static double viscous_strain_yz(const double x,
                                    const double y) {
        return 0.0;
    } // viscous_strain_yz

    static double total_strain_xx(const double x,
                                  const double y) {
        return 0.0;
    } // total_strain_xx

    static double total_strain_yy(const double x,
                                  const double y) {
        return 0.0;
    } // total_strain_yy

    static double total_strain_zz(const double x,
                                  const double y) {
        return 0.0;
    } // total_strain_zz

    static double total_strain_xy(const double x,
                                  const double y) {
        return 0.0;
    } // total_strain_xy

    static double total_strain_xz(const double x,
                                  const double y) {
        return 0.0;
    } // total_strain_xz

    static double total_strain_yz(const double x,
                                  const double y) {
        return 0.0;
    } // total_strain_yz

    static const char* strain_units(void) {
        return "none";
    } // strain_units

    // Solution subfields (nondimensional)

    // Displacement
    static double disp_x(const double x,
                         const double y,
                         const double t) {
        return 0.0;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double t) {
        return x;
    }

    // Velocity
    static double vel_x(const double x,
                        const double y,
                        const double t) {
        return 0.0;
    } // vel_x

    static double vel_y(const double x,
                        const double y,
                        const double t) {
        return 0.0;
    } // vel_y

    // Pressure
    static double fluid_pressure(const double x,
                                 const double y) {
        return 0.0;
    } // fluid_pressure

    // Trace strain
    static double trace_strain(const double x,
                               const double y,
                               const double t) {
        return 0.0;
    } // trace_strain

    static double trace_strain_dot(const double x,
                                   const double y,
                                   const double t) {
        return 0.0;
    } // trace_strain_dot

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        assert(2 == spaceDim);
        assert(x);
        assert(s);

        s[0] = disp_x(x[0], x[1], t);
        s[1] = disp_y(x[0], x[1], t);

        return 0;
    } // solnkernel_disp

    static PetscErrorCode solnkernel_fixed_disp(PetscInt spaceDim,
                                                PetscReal t,
                                                const PetscReal x[],
                                                PetscInt numComponents,
                                                PetscScalar* s,
                                                void* context) {
        assert(2 == spaceDim);
        assert(x);
        assert(s);

        s[0] = 0.0;
        s[1] = 0.0;

        return 0;
    } // solnkernel_fixed_disp

    static PetscErrorCode solnkernel_fluid_pressure(PetscInt spaceDim,
                                                    PetscReal t,
                                                    const PetscReal x[],
                                                    PetscInt numComponents,
                                                    PetscScalar* s,
                                                    void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = fluid_pressure(x[0], x[1]);

        return 0;
    } // solnkernel_fluid_pressure

    static PetscErrorCode solnkernel_trace_strain(PetscInt spaceDim,
                                                  PetscReal t,
                                                  const PetscReal x[],
                                                  PetscInt numComponents,
                                                  PetscScalar* s,
                                                  void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = trace_strain(x[0], x[1], t);

        return 0;
    } // solnkernel_trace_strain

    static PetscErrorCode solnkernel_vel(PetscInt spaceDim,
                                         PetscReal t,
                                         const PetscReal x[],
                                         PetscInt numComponents,
                                         PetscScalar* s,
                                         void* context) {
        assert(2 == spaceDim);
        assert(2 == numComponents);
        assert(s);

        s[0] = vel_x(x[0], x[1], t);
        s[1] = vel_y(x[0], x[1], t);

        return 0;
    } // solnkernel_vel

    static PetscErrorCode solnkernel_fluid_pressure_dot(PetscInt spaceDim,
                                                        PetscReal t,
                                                        const PetscReal x[],
                                                        PetscInt numComponents,
                                                        PetscScalar* s,
                                                        void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = 0.0;

        return 0;
    } // solnkernel_fluid_pressure_dot

    static PetscErrorCode solnkernel_trace_strain_dot(PetscInt spaceDim,
                                                      PetscReal t,
                                                      const PetscReal x[],
                                                      PetscInt numComponents,
                                                      PetscScalar* s,
                                                      void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = trace_strain_dot(x[0], x[1], t);

        return 0;
    } // solnkernel_trace_strain_dot

public:

    static
    TestLinearPoroviscoelasticity_Data* createData(void) {
        TestLinearPoroviscoelasticity_Data* data = new TestLinearPoroviscoelasticity_Data();assert(data);

        data->journalName = "TerzahgiPVE";
        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        // data->normalizer.setLengthScale(LENGTH_SCALE);
        // data->normalizer.setTimeScale(TIME_SCALE);
        // data->normalizer.setPressureScale(PRESSURE_SCALE);
        // data->normalizer.computeDensityScale();

        // solnDiscretizations set in derived class.

        // Material information
        data->numAuxSubfields = 12;
        static const char* _auxSubfields[12] = { // order must match order of subfields in auxiliary field
            "solid_density",
            "fluid_density",
            "fluid_viscosity",
            "porosity",
            "shear_modulus",
            "drained_bulk_modulus",
            "biot_coefficient",
            "biot_modulus",
            "maxwell_time",
            "viscous_strain",
            "total_strain",
            "isotropic_permeability",
        };
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[12] = {
            pylith::topology::Field::Discretization(1, 2), // solid_density
            pylith::topology::Field::Discretization(1, 2), // fluid_density
            pylith::topology::Field::Discretization(1, 2), // fluid_viscosity
            pylith::topology::Field::Discretization(1, 2), // porosity
            pylith::topology::Field::Discretization(1, 2), // shear_modulus
            pylith::topology::Field::Discretization(1, 2), // drained_bulk_modulus
            pylith::topology::Field::Discretization(1, 2), // biot_coefficient
            pylith::topology::Field::Discretization(1, 2), // biot_modulus
            pylith::topology::Field::Discretization(1, 2), // maxwell_time
            pylith::topology::Field::Discretization(1, 2), // viscous_strain
            pylith::topology::Field::Discretization(1, 2), // total_strain
            pylith::topology::Field::Discretization(1, 2), // isotropic_permeability
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        data->auxDB.addValue("solid_density", solid_density, density_units());
        data->auxDB.addValue("fluid_density", fluid_density, density_units());
        data->auxDB.addValue("fluid_viscosity", fluid_viscosity, viscosity_units());
        data->auxDB.addValue("porosity", porosity, porosity_units());
        data->auxDB.addValue("shear_modulus", shear_modulus, modulus_units());
        data->auxDB.addValue("drained_bulk_modulus", drained_bulk_modulus, modulus_units());
        data->auxDB.addValue("biot_coefficient", biot_coefficient, modulus_units());
        data->auxDB.addValue("fluid_bulk_modulus", fluid_bulk_modulus, modulus_units());
        data->auxDB.addValue("solid_bulk_modulus", solid_bulk_modulus, modulus_units());
        data->auxDB.addValue("solid_viscosity", solid_viscosity, viscosity_units());

        data->auxDB.addValue("viscous_strain_xx", viscous_strain_xx, strain_units());
        data->auxDB.addValue("viscous_strain_yy", viscous_strain_yy, strain_units());
        data->auxDB.addValue("viscous_strain_zz", viscous_strain_zz, strain_units());
        data->auxDB.addValue("viscous_strain_xy", viscous_strain_xy, strain_units());
        data->auxDB.addValue("viscous_strain_xz", viscous_strain_xz, strain_units());
        data->auxDB.addValue("viscous_strain_yz", viscous_strain_yz, strain_units());

        data->auxDB.addValue("total_strain_xx", total_strain_xx, strain_units());
        data->auxDB.addValue("total_strain_yy", total_strain_yy, strain_units());
        data->auxDB.addValue("total_strain_zz", total_strain_zz, strain_units());
        data->auxDB.addValue("total_strain_xy", total_strain_xy, strain_units());
        data->auxDB.addValue("total_strain_xz", total_strain_xz, strain_units());
        data->auxDB.addValue("total_strain_yz", total_strain_yz, strain_units());

        data->auxDB.addValue("isotropic_permeability", isotropic_permeability, permeability_units());
        data->auxDB.setCoordSys(data->cs);

        data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
        data->rheology.useReferenceState(false);

        data->material.setIdentifier("poroviscoelasticity");
        data->material.setName("material-id=24");
        data->material.setLabelValue(24);

        static const PylithInt constrainedXY[2] = { 0, 1 };
        static const PylithInt numConstrained = 2;
        data->bcs.resize(4);
        { // Displacement -x
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedXY, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[0] = bc;
        }
        { // Displacement +x
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedXY, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[1] = bc;
        }
        { // Displacement -ys
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_yneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedXY, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[2] = bc;
        }
        { // Displacement +y
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_ypos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedXY, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[3] = bc;
        }

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[3] = {
            solnkernel_disp,
            solnkernel_fluid_pressure,
            solnkernel_trace_strain,
        };

        static const pylith::testing::MMSTest::solution_fn _exactSolnDotFns[3] = {
            solnkernel_vel,
            solnkernel_fluid_pressure_dot,
            solnkernel_trace_strain_dot,
        };

        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        data->exactSolnDotFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnDotFns);

        return data;

    } // createData

}; // SimpleShearPVE
const double pylith::_SimpleShearPVE::LENGTH_SCALE = 1.0e+3;
const double pylith::_SimpleShearPVE::TIME_SCALE = 2.0;
const double pylith::_SimpleShearPVE::PRESSURE_SCALE = 2.25e+10;

const double pylith::_SimpleShearPVE::PRESSURE0 = 4.0e+6;
const double pylith::_SimpleShearPVE::XMAX = 8.0e+3;

// ------------------------------------------------------------------------------------------------
pylith::TestLinearPoroviscoelasticity_Data*
pylith::SimpleShearPVE::QuadQ2Q1Q1(void) {
    TestLinearPoroviscoelasticity_Data* data = pylith::_SimpleShearPVE::createData();assert(data);

    data->meshFilename = "data/quad_simple.mesh";

    data->numSolnSubfields = 3;
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(1, 2), // fluid pressure
        pylith::topology::Field::Discretization(1, 2), // trace strain
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[12] = {
        pylith::topology::Field::Discretization(1, 2), // solid_density
        pylith::topology::Field::Discretization(1, 2), // fluid_density
        pylith::topology::Field::Discretization(1, 2), // fluid_viscosity
        pylith::topology::Field::Discretization(1, 2), // porosity
        pylith::topology::Field::Discretization(1, 2), // shear_modulus
        pylith::topology::Field::Discretization(1, 2), // drained_bulk_modulus
        pylith::topology::Field::Discretization(1, 2), // biot_coefficient
        pylith::topology::Field::Discretization(1, 2), // biot_modulus
        pylith::topology::Field::Discretization(1, 2), // maxwell_time
        pylith::topology::Field::Discretization(1, 2), // viscous_strain
        pylith::topology::Field::Discretization(1, 2), // total_strain
        pylith::topology::Field::Discretization(1, 2), // isotropic_permeability
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ2Q1Q1


// End of file
