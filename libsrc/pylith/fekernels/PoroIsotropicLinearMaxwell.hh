/*
 * ================================================================================================
 * This code is part of PyLith, developed through the Computational Infrastructure
 * for Geodynamics (https://github.com/geodynamics/pylith).
 *
 * Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
 * All rights reserved.
 *
 * See https://mit-license.org/ and LICENSE.md and for license information.
 * =================================================================================================
 */
#pragma once

/*
 * Kernels for linear poroelasticity plane strain.
 *
 * Solution fields: [disp(dim), pressure(1),trace_strain(1) ] (QS)
 * OR
 * Solution fields: [disp(dim), pressure(1),velocity(dim) ] (DYN)
 *
 * Auxiliary fields:
 * -- numA : number of auxiliary fields
 ***** Required fields(govening equations) + option fields + required fields (rheology)
 * - 0: solid_density(1)
 * - 1: fluid_density(1)
 * - 2: fluid_viscosity(1)
 * - 3: porosity(1)
 *
 ** Optional fields
 * - +1: gravity_field (dim, optional) (4)
 * - +1: body_force(dim,optional) (4,5)
 * - +1: source_density(1,optional) (4,5,6)
 * - +1: reference_stress(4,optional) (stress_xx, stress_yy, stress_xy, stress_zz)  // numA - 7
 *     2D: 4 components (stress_xx, stress_yy, stress_zz, stress_xy)
 *     3D: 6 components (stress_xx, stress_yy, stress_zz, stress_xy, stress_yz, stress_xz)
 * - +1: reference_strain(4,optional) (strain_xx, strain_yy, strain_xy, strain_zz) // numA - 6
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 *
 ** Rheological fields
 * - numA - 5: addShearModulus(1)
 * - numA - 4: addDrainedBulkModulus(1)
 * - numA - 3: addBiotCoefficient(1)
 *      Isotropic: numA - 2: addIsotropicPermeability(1)
 *      Tensor: numA - 2: addTensorPermeability(4,optional) (permeability_xx, permeability_yy, permeability_zz,
 * permeability_xy)
 *          2D: 4 components (permeability_xx, permeability_yy, permeability_zz, permeability_xy)
 *          3D: 6 components (permeability_xx, permeability_yy, permeability_zz, permeability_xy, permeability_yz,
 * permeability_xz)
 * - numA - 1: addFluidBulkModulus(1)
 *
 * The poroelasticity subfields come first (with required ones before optional ones) followed by the rheology subfields
 * (optional ones before required ones). The rheology fields have required fields last because we index from the back.
 *
 * :TODO: @robert Add equation here
 *
 *
 * ======================================================================
 *
 * Kernel interface.
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field.
 * @param[in] numA Number of registered subfields in auxiliary field.
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] numConstants Number of registered constants.
 * @param[in] constants Array of registered constants.
 * @param[out] f0 [dim].
 */

// Include directives ---------------------------------------------------
#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations
#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/IsotropicLinearMaxwell.hh" // USES Isotropic Linear Maxwell context and kernels
#include "pylith/fekernels/IsotropicLinearPoroelasticity.hh" //Uses Isotropic Linear Poroelasticity kernels 

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linear poroelasticity.
class pylith::fekernels::PoroIsotropicLinearMaxwell {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

   // --------------------------------------------------------------------------------------------
    static inline
    void setIsotropicLinearPoroelasticityContext(pylith::fekernels::IsotropicLinearPoroelasticity::Context* context,
                                                const PylithInt dim,
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
                                                const PylithScalar x[],
                                                const PylithInt numConstants,
                                                const PylithScalar constants[],
                                                const pylith::fekernels::TensorOps& tensorOps) {
            assert(context);

            // Incoming solution subfields
            const PylithInt i_pressure = 1;

            // Incoming poroelastic auxiliary subfields
            const PylithInt i_fluidViscosity = 2;

            // Incoming rheology auxiliary subfields.
            const PylithInt i_shearModulus = numA - 8;
            const PylithInt i_drainedBulkModulus = numA - 7;
            const PylithInt i_biotCoefficient = numA - 6;
            const PylithInt i_biotModulus = numA - 5;

            assert(numA >= 11); // also have density
            assert(s);
            assert(sOff);
            assert(sOff[i_pressure] >= 0);
            assert(a);
            assert(aOff);
            assert(aOff[i_shearModulus] >= 0);
            assert(aOff[i_drainedBulkModulus] >= 0);
            assert(aOff[i_biotCoefficient] >= 0);
            assert(aOff[i_biotModulus] >= 0);
            assert(aOff[i_fluidViscosity] >= 0);
            assert(constants);

             // Solution Variables
            context->pressure = s[sOff[i_pressure]];

            // Poroelastic Auxiliary Variables
            context->fluidViscosity = a[aOff[i_fluidViscosity]];assert(context->fluidViscosity > 0.0);

            // Rheology Specific Auxiliary Variables
            context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
            context->drainedBulkModulus = a[aOff[i_drainedBulkModulus]];assert(context->drainedBulkModulus > 0.0);
            context->biotCoefficient = a[aOff[i_biotCoefficient]];assert(context->biotCoefficient > 0.0);
            context->biotModulus = a[aOff[i_biotModulus]];assert(context->biotModulus > 0.0);
        
       

    } // setIsotropicLinearPoroelasticityContext

    static inline
    void setIsotropicLinearMaxwellContext(pylith::fekernels::IsotropicLinearMaxwell::Context* context,
                                         const PylithInt dim,
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
                                         const PylithScalar x[],
                                         const PylithInt numConstants,
                                         const PylithScalar constants[],
                                         const pylith::fekernels::TensorOps& tensorOps) {
            assert(context);

            // Incoming rheology auxiliary subfields.
            const PylithInt i_shearModulus = numA - 8;
            const PylithInt i_drainedBulkModulus = numA - 7;
            const PylithInt i_maxwellTime = numA - 4;
            const PylithInt i_viscousStrain = numA - 3;
            const PylithInt i_totalStrain = numA - 2;

            assert(numA >= 11); // also have density
            assert(s);
            assert(sOff);
            assert(aOff[i_maxwellTime] >= 0);
            assert(aOff[i_viscousStrain] >= 0);
            assert(aOff[i_totalStrain] >= 0);
            assert(1 == numConstants);
            assert(constants);

            // Rheology Specific Auxiliary Variables
            context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
            context->bulkModulus = a[aOff[i_drainedBulkModulus]];assert(context->bulkModulus > 0.0); //set IsotropicLinearMaxwell context to have its bulk modulus be the drained bulk modulus
            context->maxwellTime = a[aOff[i_maxwellTime]];assert(context->maxwellTime > 0.0);
            context->dt = constants[0];

            tensorOps.fromVector(&a[aOff[i_viscousStrain]], &context->viscousStrain);
            tensorOps.fromVector(&a[aOff[i_totalStrain]], &context->totalStrain);
        } //setIsotropicLinearMaxwellContext

     // --------------------------------------------------------------------------------------------
    static inline
    void setIsotropicLinearPoroelasticityContextIsotropicPerm(pylith::fekernels::IsotropicLinearPoroelasticity::Context* context,
                                                                const PylithInt dim,
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
                                                                const PylithScalar x[],
                                                                const PylithInt numConstants,
                                                                const PylithScalar constants[],
                                                                const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);

        // Incoming auxiliary fields.
        const PylithInt i_isotropicPermeability = numA - 1;

        // Using isotropic permeability
        tensorOps.fromScalar(a[aOff[i_isotropicPermeability]], &context->permeability);

    } // setIsotropicLinearPoroelasticityContextIsotropicPerm

       // --------------------------------------------------------------------------------------------
    static inline
    void setIsotropicLinearPoroelasticityContextTensorPerm(pylith::fekernels::IsotropicLinearPoroelasticity::Context* context,
                                                            const PylithInt dim,
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
                                                            const PylithScalar x[],
                                                            const PylithInt numConstants,
                                                            const PylithScalar constants[],
                                                            const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);

        // Incoming auxiliary fields.
        const PylithInt i_tensorPermeability = numA - 1;

        // Using tensor permeability
        tensorOps.fromVector(&a[aOff[i_tensorPermeability]], &context->permeability);

    } // setIsotropicLinearPoroelasticityContextTensorPerm

        // --------------------------------------------------------------------------------------------
    static inline
    void setIsotropicLinearPoroelasticityContextRefState(pylith::fekernels::IsotropicLinearPoroelasticity::Context* context,
                                                    const PylithInt dim,
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
                                                    const PylithScalar x[],
                                                    const PylithInt numConstants,
                                                    const PylithScalar constants[],
                                                    const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);
        
        // Incoming auxiliary fields.
        const PylithInt i_refStress = numA - 10;
        const PylithInt i_refStrain = numA - 9;

        // Reference stress and strain
        tensorOps.fromVector(&a[aOff[i_refStress]], &context->refStress);
        tensorOps.fromVector(&a[aOff[i_refStrain]], &context->refStrain);

    } // setIsotropicLinearPoroelasticityContextRefState

    // --------------------------------------------------------------------------------------------
    static inline
    void setIsotropicLinearMaxwellContextRefState(pylith::fekernels::IsotropicLinearMaxwell::Context* context,
                                                        const PylithInt dim,
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
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);
        
        // Incoming auxiliary fields.
        const PylithInt i_refStress = numA - 10;
        const PylithInt i_refStrain = numA - 9;
        const PylithInt i_shearModulus = numA - 8;
        const PylithInt i_drainedBulkModulus = numA - 7;
        const PylithInt i_maxwellTime = numA - 4;
        const PylithInt i_viscousStrain = numA - 3;
        const PylithInt i_totalStrain = numA - 2;

        assert(numA >= 11); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_drainedBulkModulus] >= 0);
        assert(aOff[i_maxwellTime] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_totalStrain] >= 0);
        assert(1 == numConstants);
        assert(constants);

        context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
        context->bulkModulus = a[aOff[i_drainedBulkModulus]];assert(context->bulkModulus > 0.0);
        context->maxwellTime = a[aOff[i_maxwellTime]];assert(context->maxwellTime > 0.0);
        context->dt = constants[0];

        tensorOps.fromVector(&a[aOff[i_viscousStrain]], &context->viscousStrain);
        tensorOps.fromVector(&a[aOff[i_totalStrain]], &context->totalStrain);
        tensorOps.fromVector(&a[aOff[i_refStress]], &context->refStress);
        tensorOps.fromVector(&a[aOff[i_refStrain]], &context->refStrain);

    } // setIsotropicLinearMaxwellContextRefState


    // --------------------------------------------------------------------------------------------
    static inline
    void setIsotropicLinearPoroelasticityContextQuasistatic(pylith::fekernels::IsotropicLinearPoroelasticity::Context* context,
                                                            const PylithInt dim,
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
                                                            const PylithScalar x[],
                                                            const PylithInt numConstants,
                                                            const PylithScalar constants[],
                                                            const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);

        // Incoming solution fields.
        const PylithInt i_trace_strain = 2;
        assert(sOff[i_trace_strain] >= 0);

        // Variables &c
        context->trace_strain = s[sOff[i_trace_strain]];

    } // setIsotropicLinearPoroelasticityContextQuasistatic
   

    // --------------------------------------------------------------------------------------------
    /** Helper function for calculating Cauchy stress for WITHOUT a reference stress and strain.
     *
     * ISA Poroelasticity::stressFn
     *
     * @param[in] isotropicLinearPoroelasticityContext IsotropicLinearPoroelasticity context.
     * @param[in] isotropicLinearMaxwellContext IsotropicLinearMaxwell context.
     * @param[in] strain Strain tensor.
     * @param[in] tensorOps Tensor operations.
     * @param[out] stress Stress tensor.
     */
    static inline
    void cauchyStress(const pylith::fekernels::IsotropicLinearPoroelasticity::Context* isotropicLinearPoroelasticityContext,
                      const pylith::fekernels::IsotropicLinearMaxwell::Context* isotropicLinearMaxwellContext,
                      const pylith::fekernels::Tensor& strain,
                      const pylith::fekernels::TensorOps& tensorOps,
                      pylith::fekernels::Tensor* stress) {
        assert(isotropicLinearPoroelasticityContext);
        assert(isotropicLinearMaxwellContext);
        assert(stress);

        pylith::fekernels::IsotropicLinearPoroelasticity::meanStress(isotropicLinearPoroelasticityContext->pressure, isotropicLinearPoroelasticityContext->trace_strain, isotropicLinearPoroelasticityContext->drainedBulkModulus, 
                                                                     isotropicLinearPoroelasticityContext->biotCoefficient, strain, stress);

        const PylithReal dt = isotropicLinearMaxwellContext->dt;assert(dt > 0.0);
        const PylithReal maxwellTime = isotropicLinearMaxwellContext->maxwellTime;
        const pylith::fekernels::Tensor& totalStrain = isotropicLinearMaxwellContext->totalStrain;
        const pylith::fekernels::Tensor& viscousStrainPrev = isotropicLinearMaxwellContext->viscousStrain;
        pylith::fekernels::Tensor viscousStrain;
        pylith::fekernels::IsotropicLinearMaxwell::viscousStrain(maxwellTime, viscousStrainPrev, totalStrain, strain, dt, &viscousStrain);

        deviatoricStress(isotropicLinearPoroelasticityContext->trace_strain, isotropicLinearMaxwellContext->shearModulus, viscousStrain, stress);

    } // cauchyStress

    // --------------------------------------------------------------------------------------------
    /** Helper function for calculating Cauchy stress WITH reference stress/strain.
     *
     * @param[in] isotropicLinearPoroelasticityContext IsotropicLinearPoroelasticity context.
     * @param[in] isotropicLinearMaxwellContext IsotropicLinearMaxwell context.
     * @param[in] strain Strain tensor.
     * @param[in] tensorOps Tensor operations.
     * @param[out] stress Stress tensor.
     */
    static inline
    void cauchyStress_refState(pylith::fekernels::IsotropicLinearPoroelasticity::Context* isotropicLinearPoroelasticityContext,
                      pylith::fekernels::IsotropicLinearMaxwell::Context* isotropicLinearMaxwellContext,
                      const pylith::fekernels::Tensor& strain,
                      const pylith::fekernels::TensorOps& tensorOps,
                      pylith::fekernels::Tensor* stress) {
        assert(isotropicLinearPoroelasticityContext);
        assert(isotropicLinearMaxwellContext);
        assert(stress);

        const pylith::fekernels::Tensor& refStress = isotropicLinearPoroelasticityContext->refStress;
        const pylith::fekernels::Tensor& refStrain = isotropicLinearPoroelasticityContext->refStrain;
        pylith::fekernels::IsotropicLinearPoroelasticity::meanStress_refState(isotropicLinearPoroelasticityContext->pressure, isotropicLinearPoroelasticityContext->trace_strain, isotropicLinearPoroelasticityContext->drainedBulkModulus,
                                                                              isotropicLinearPoroelasticityContext->biotCoefficient, refStress, refStrain, strain, stress);

        const PylithReal dt = isotropicLinearMaxwellContext->dt;assert(dt > 0.0);
        const PylithReal maxwellTime = isotropicLinearMaxwellContext->maxwellTime;
        const pylith::fekernels::Tensor& totalStrain = isotropicLinearMaxwellContext->totalStrain;
        const pylith::fekernels::Tensor& viscousStrainPrev = isotropicLinearMaxwellContext->viscousStrain;
        pylith::fekernels::Tensor viscousStrain;
        pylith::fekernels::IsotropicLinearMaxwell::viscousStrain(maxwellTime, viscousStrainPrev, totalStrain, strain, dt, &viscousStrain);

        deviatoricStress_refState(isotropicLinearPoroelasticityContext->trace_strain, isotropicLinearMaxwellContext->shearModulus, refStress, refStrain, viscousStrain, stress);
    } // cauchyStress_refState

     // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITHOUT reference stress and strain.
     */
    static inline
    void deviatoricStress(const PylithReal strainTrace,
                          const PylithReal shearModulus,
                          const pylith::fekernels::Tensor& viscousStrain,
                          pylith::fekernels::Tensor* stress) {
        assert(shearModulus > 0.0);
        assert(stress);

        // const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        const PylithReal traceTerm = -2.0/3.0*shearModulus * strainTrace;

        stress->xx += 2.0*shearModulus*viscousStrain.xx + traceTerm;
        stress->yy += 2.0*shearModulus*viscousStrain.yy + traceTerm;
        stress->zz += 2.0*shearModulus*viscousStrain.zz + traceTerm;
        stress->xy += 2.0*shearModulus*viscousStrain.xy;
        stress->yz += 2.0*shearModulus*viscousStrain.yz;
        stress->xz += 2.0*shearModulus*viscousStrain.xz;
    } // deviatoricStress

    static inline
    void deviatoricStress_refState(const PylithReal strainTrace,
                                   const PylithReal shearModulus,
                                   const pylith::fekernels::Tensor& refStress,
                                   const pylith::fekernels::Tensor& refStrain,
                                   const pylith::fekernels::Tensor& viscousStrain,
                                   pylith::fekernels::Tensor* stress) {
        assert(shearModulus > 0.0);
        assert(stress);

        // const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        const PylithReal refStrainTrace = refStrain.xx + refStrain.yy + refStrain.zz;
        const PylithReal meanRefStress = (refStress.xx + refStress.yy + refStress.zz) / 3.0;
        const PylithReal traceTerm = -2.0/3.0*shearModulus * (strainTrace - refStrainTrace);

        stress->xx += refStress.xx - meanRefStress + 2.0*shearModulus*(viscousStrain.xx-refStrain.xx) + traceTerm;
        stress->yy += refStress.yy - meanRefStress + 2.0*shearModulus*(viscousStrain.yy-refStrain.yy) + traceTerm;
        stress->zz += refStress.zz - meanRefStress + 2.0*shearModulus*(viscousStrain.zz-refStrain.zz) + traceTerm;
        stress->xy += refStress.xy + 2.0*shearModulus*(viscousStrain.xy - refStrain.xy);
        stress->yz += refStress.yz + 2.0*shearModulus*(viscousStrain.yz - refStrain.yz);
        stress->xz += refStress.xz + 2.0*shearModulus*(viscousStrain.xz - refStrain.xz);
    }

}; //PoroIsotropicLinearMaxwell

/// Kernels specific to isotropic, linear poroelasticity in Plane Strain.
class pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain {
        // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    // ================================= LHS =======================================

     // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms.
    static inline
    void f0p_implicit(const PylithInt dim,
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
                      const PylithScalar x[],
                      const PylithInt numConstants,
                      const PylithScalar constants[],
                      PylithScalar f0[]) {
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit_context(
            _dim, s_t, &poroelasticContext, &rheologyContext, f0);

    } // f0p_implicit

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source(const PylithInt dim,
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
                             const PylithScalar x[],
                             const PylithInt numConstants,
                             const PylithScalar constants[],
                             PylithScalar f0[]) {
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit_source_context(
            _dim, &poroelasticContext, &rheologyContext, f0);

    } // f0p_implicit_source

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static
    void f0p_implicit_source_body(const PylithInt dim,
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
                                  const PylithScalar x[],
                                  const PylithInt numConstants,
                                  const PylithScalar constants[],
                                  PylithScalar f0[]) {
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextBodyForceSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit_source_body_context(
            _dim, &poroelasticContext, &rheologyContext, f0);
    } // f0p_implicit_source_body

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_grav(const PylithInt dim,
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
                                  const PylithScalar x[],
                                  const PylithInt numConstants,
                                  const PylithScalar constants[],
                                  PylithScalar f0[]) {
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextGravitySourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit_source_grav_context(
            _dim, &poroelasticContext, &rheologyContext, f0);
       
    } // f0p_implicit_source_grav

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_grav_body(const PylithInt dim,
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
                                       const PylithScalar x[],
                                       const PylithInt numConstants,
                                       const PylithScalar constants[],
                                       PylithScalar f0[]) {
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextGravityBodyForceSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit_source_grav_body_context(
            _dim, &poroelasticContext, &rheologyContext, f0);
    } // f0p_implicit_source_grav_body

     // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear visco-poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Quasi - Static Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static inline
    void f1u(const PylithInt dim,
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
             const PylithScalar x[],
             const PylithInt numConstants,
             const PylithScalar constants[],
             PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContextIsotropicLinearPoroelasticity;
        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContextIsotropicLinearMaxwell;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextQuasistatic(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContext(
            &rheologyContextIsotropicLinearMaxwell, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // pylith::fekernels::Elasticity::f1v(
        //     strainContext, &rheologyContextIsotropicLinearPoroelasticity,
        //     pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
        //     pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress,
        //     pylith::fekernels::Tensor::ops2D,
        //     f1);

        Tensor strain;
        pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain(strainContext, &strain);

        Tensor stress;
        TensorOps tensorOps = pylith::fekernels::Tensor::ops2D;
        pylith::fekernels::PoroIsotropicLinearMaxwell::cauchyStress(&rheologyContextIsotropicLinearPoroelasticity, &rheologyContextIsotropicLinearMaxwell, strain, tensorOps, &stress);

        PylithScalar stressTensor[9] = {0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
        tensorOps.toTensor(stress, stressTensor);

        for (PylithInt i = 0; i < _dim*_dim; ++i) {
                f1[i] -= stressTensor[i];
        } // for

    } // f1u

     // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static inline
    void f1u_refstate(const PylithInt dim,
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
                      const PylithScalar x[],
                      const PylithInt numConstants,
                      const PylithScalar constants[],
                      PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContextIsotropicLinearPoroelasticity;
        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContextIsotropicLinearMaxwell;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextQuasistatic(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextRefState(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContext(
            &rheologyContextIsotropicLinearMaxwell, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContextRefState(
            &rheologyContextIsotropicLinearMaxwell, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        Tensor strain;
        pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain(strainContext, &strain);

        Tensor stress;
        TensorOps tensorOps = pylith::fekernels::Tensor::ops2D;
        pylith::fekernels::PoroIsotropicLinearMaxwell::cauchyStress_refState(&rheologyContextIsotropicLinearPoroelasticity, &rheologyContextIsotropicLinearMaxwell, strain, tensorOps, &stress);

        PylithScalar stressTensor[9] = {0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
        tensorOps.toTensor(stress, stressTensor);

        for (PylithInt i = 0; i < _dim*_dim; ++i) {
            f1[i] -= stressTensor[i];
        } // for


    } // f1u_refstate

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity
     *
     */
    static inline
    void f1p(const PylithInt dim,
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
             const PylithScalar x[],
             const PylithInt numConstants,
             const PylithScalar constants[],
             PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_tensor_permeability(const PylithInt dim,
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
                                 const PylithScalar x[],
                                 const PylithInt numConstants,
                                 const PylithScalar constants[],
                                 PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body(const PylithInt dim,
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
                  const PylithScalar x[],
                  const PylithInt numConstants,
                  const PylithScalar constants[],
                  PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_body

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_tensor_permeability(const PylithInt dim,
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
                                      const PylithScalar x[],
                                      const PylithInt numConstants,
                                      const PylithScalar constants[],
                                      PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_body_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_gravity(const PylithInt dim,
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
                     const PylithScalar x[],
                     const PylithInt numConstants,
                     const PylithScalar constants[],
                     PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_gravity

     // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_gravity_tensor_permeability(const PylithInt dim,
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
                                         const PylithScalar x[],
                                         const PylithInt numConstants,
                                         const PylithScalar constants[],
                                         PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_gravity_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_gravity(const PylithInt dim,
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
                          const PylithScalar x[],
                          const PylithInt numConstants,
                          const PylithScalar constants[],
                          PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravityBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_body_gravity

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_gravity_tensor_permeability(const PylithInt dim,
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
                                              const PylithScalar x[],
                                              const PylithInt numConstants,
                                              const PylithScalar constants[],
                                              PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravityBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_body_gravity_tensor_permeability

    // =========================== LHS Jacobian ============================

    // ----------------------------------------------------------------------
    /* Jf3_uu entry function for isotropic linear visco-poroelasticity.
     */
    static inline
    void Jf3uu(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf3[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context isotropicLinearPoroelasticityRheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &isotropicLinearPoroelasticityRheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearMaxwell::Context isotropicLinearMaxwellRheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContext(
            &isotropicLinearMaxwellRheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        const PylithScalar shearModulus = isotropicLinearPoroelasticityRheologyContext.shearModulus;
        const PylithScalar maxwellTime = isotropicLinearMaxwellRheologyContext.maxwellTime;
        const PylithScalar dt = isotropicLinearMaxwellRheologyContext.dt;

        const PylithScalar dq = pylith::fekernels::IsotropicLinearMaxwell::viscousStrainCoeff(dt, maxwellTime);

        //Unique components of Jacobian.
        const PylithReal C1111 = 4.0/3.0 * shearModulus * dq;
        const PylithReal C1122 = -2.0/3.0 * shearModulus * dq;
        const PylithReal C1212 = shearModulus * dq;

        /* j(f,g,df,dg) = C(f,df,g,dg)
         *
         * 0:  j0000 = C1111 = 1.0*drainedBulkModulus + 1.33333333333333*delHM*shearModulus
         * 1:  j0001 = C1112 = 0
         * 2:  j0010 = C1211 = 0
         * 3:  j0011 = C1212 = 1.0*delHM*shearModulus
         * 4:  j0100 = C1121 = 0
         * 5:  j0101 = C1122 = 1.0*drainedBulkModulus - 0.666666666666667*delHM*shearModulus
         * 6:  j0110 = C1221 = 1.0*delHM*shearModulus
         * 7:  j0111 = C1222 = 0
         * 8:  j1000 = C2111 = 0
         * 9:  j1001 = C2112 = 1.0*delHM*shearModulus
         * 10:  j1010 = C2211 = 1.0*drainedBulkModulus - 0.666666666666667*delHM*shearModulus
         * 11:  j1011 = C2212 = 0
         * 12:  j1100 = C2121 = 1.0*delHM*shearModulus
         * 13:  j1101 = C2122 = 0
         * 14:  j1110 = C2221 = 0
         * 15:  j1111 = C2222 = 1.0*drainedBulkModulus + 1.33333333333333*delHM*shearModulus
         */

        /* Nonzero Jacobian entries. */
        Jf3[0] -= C1111;// - 2.6666 * shearModulus; /* j0000 */
        Jf3[3] -= C1212;// - shearModulus; /* j0011 */
        Jf3[5] -= C1122; /* j0101 */
        Jf3[6] -= C1212;// - shearModulus; /* j0110 */
        Jf3[9] -= C1212;// - shearModulus; /* j1001 */
        Jf3[10] -= C1122; /* j1010 */
        Jf3[12] -= C1212;// - shearModulus; /* j1100 */
        Jf3[15] -= C1111;// - 2.6666 * shearModulus; /* j1111 */

    }

    // ----------------------------------------------------------------------
    /** Jf2_up entry function for isotropic linear poroelasticity.
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf2up(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf2[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf2up_context(
            dim, &rheologyContext, Jf2);
    } // Jf2up

    // -----------------------------------------------------------------------------
    // Jf2ue function for isotropic linear poroviscoelasticity.
    static inline
    void Jf2ue(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf2[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf2ue_context(
            dim, &rheologyContext, Jf2);
    } // Jf2ue

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroviscoelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf3pp(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf3[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf3pp_context(
            dim, &rheologyContext, Jf3);

    } // Jf3pp

       // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroviscoelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf3pp_tensor_permeability(const PylithInt dim,
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
                                   const PylithReal s_tshift,
                                   const PylithScalar x[],
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar Jf3[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf3pp_context(
            dim, &rheologyContext, Jf3);

    } // Jf3pp_tensor_permeability

    // ----------------------------------------------------------------------
    /** Jf0_pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0pp(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf0[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf0pp_context(
            _dim, s_tshift, &rheologyContext, Jf0);
        
    } // Jf0pp

     // ----------------------------------------------------------------------
    /** Jf0_pe entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0pe(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf0[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf0pe_context(
            _dim, s_tshift, &rheologyContext, Jf0);

    } // Jf0pe

     // ===========================================================================================
    // Kernels for updating state variables
    // ===========================================================================================

    // Use pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector() to update total strain.

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating viscous strain as a vector for 2D plane strain isotropic
     * linear Maxwell viscoelasticity.
     *
     * Used to output viscous strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void viscousStrain_infinitesimalStrain_asVector(const PylithInt dim,
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
                                                    const PylithScalar x[],
                                                    const PylithInt numConstants,
                                                    const PylithScalar constants[],
                                                    PylithScalar viscousStrain[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearMaxwell::viscousStrain_asVector(
            strainContext, rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::Tensor::ops2D,
            viscousStrain);
    }

// ---------------------------------------------------------------------------------------------------------------------
    /* Update porosity for a linear poroelastic material, implicit.
     */
    static inline
    void updatePorosityImplicit(const PylithInt dim,
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
                                const PylithScalar x[],
                                const PylithInt numConstants,
                                const PylithScalar constants[],
                                PylithScalar porosity[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Incoming solution fields.
        const PylithInt i_pressure_t = 4;
        const PylithInt i_trace_strain_t = 5;

        // Incoming re-packed auxiliary field.

        // Poroelasticity
        const PylithInt i_porosityPrev = 3;

        // IsotropicLinearPoroelasticity
        const PylithInt i_drainedBulkModulus = numA - 7;
        const PylithInt i_biotCoefficient = numA - 6;

        // Constants
        const PylithScalar dt = constants[0];

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 6);
        assert(numA >= 6);
        assert(aOff);
        assert(aOff[i_porosityPrev] >= 0);
        assert(porosity);

        // Do stuff
        const PylithScalar pressure_t = s ? s[sOff[i_pressure_t]] : 0.0;
        const PylithScalar trace_strain_t = s ? s[sOff[i_trace_strain_t]] : 0.0;

        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar porosityPrev = a[aOff[i_porosityPrev]];

        // Update porosity
        porosity[0] = porosityPrev + dt * ((biotCoefficient - porosityPrev) * trace_strain_t +
                                           ((1.0 - biotCoefficient) * (biotCoefficient - porosityPrev)) /
                                           drainedBulkModulus * pressure_t);
        porosity[0] = std::max(0.0, std::min(1.0, porosity[0]));


    } // updatePorosityImplicit

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 2D isotropic linear poroelasticity with
     * infinitesimal strain WITHOUT reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., biot_coefficient(1), shear_modulus(1), drained_bulk_modulus(1)]
     */
    static inline
    void cauchyStress_infinitesimalStrain_asVector(const PylithInt dim,
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
                                                   const PylithScalar x[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar stressVector[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContextIsotropicLinearPoroelasticity;
        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContextIsotropicLinearMaxwell;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextQuasistatic(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContext(
            &rheologyContextIsotropicLinearMaxwell, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);


        Tensor strain;
        pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain(strainContext, &strain);

        Tensor stress;
        TensorOps tensorOps = pylith::fekernels::Tensor::ops2D;
        pylith::fekernels::PoroIsotropicLinearMaxwell::cauchyStress(&rheologyContextIsotropicLinearPoroelasticity, &rheologyContextIsotropicLinearMaxwell, strain, tensorOps, &stress);

        tensorOps.toVector(stress, stressVector);


        // pylith::fekernels::Elasticity::stress_asVector(
        //     strainContext, &rheologyContext,
        //     pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
        //     pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress,
        //     pylith::fekernels::Tensor::ops2D,
        //     stressVector);

    } // cauchyStress_infinitesimalStrain_asVector

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for DD isotropic linear poroelasticity with
     * infinitesimal strain WITH reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., biot_coefficient(1), shear_modulus(1), drained_bulk_modulus(1)]
     */
    static inline
    void cauchyStress_infinitesimalStrain_refState_asVector(const PylithInt dim,
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
                                                   const PylithScalar x[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar stressVector[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContextIsotropicLinearPoroelasticity;
        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContextIsotropicLinearMaxwell;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextRefState(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        // Using dynamic formulation for trace strain, assuming that it will be equal to the variable
        // for QS
        // pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextDynamic(
        //     &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
        //     t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContext(
            &rheologyContextIsotropicLinearMaxwell, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContextRefState(
            &rheologyContextIsotropicLinearMaxwell, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        Tensor strain;
        pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain(strainContext, &strain);

        Tensor stress;
        TensorOps tensorOps = pylith::fekernels::Tensor::ops2D;
        pylith::fekernels::PoroIsotropicLinearMaxwell::cauchyStress_refState(&rheologyContextIsotropicLinearPoroelasticity, &rheologyContextIsotropicLinearMaxwell, strain, tensorOps, &stress);

        tensorOps.toVector(stress, stressVector);

    } // cauchyStress_infinitesimalStrain_refState_asVector

        // Calculate water content
    static inline
    void waterContent_asScalar(const PylithInt dim,
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
                               const PylithScalar x[],
                               const PylithInt numConstants,
                               const PylithScalar constants[],
                               PylithReal* waterContent) {
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::waterContent_asScalar_context(
            dim, &poroelasticContext, &rheologyContext, waterContent);

    } // waterContent_asScalar

};// PoroIsotropicLinearMaxwellPlaneStrain

class pylith::fekernels::PoroIsotropicLinearMaxwell3D {
            // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    // ================================= LHS =======================================

     // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms.
    static inline
    void f0p_implicit(const PylithInt dim,
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
                      const PylithScalar x[],
                      const PylithInt numConstants,
                      const PylithScalar constants[],
                      PylithScalar f0[]) {
        const PylithInt _dim = 3;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit_context(
            _dim, s_t, &poroelasticContext, &rheologyContext, f0);

    } // f0p_implicit

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source(const PylithInt dim,
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
                             const PylithScalar x[],
                             const PylithInt numConstants,
                             const PylithScalar constants[],
                             PylithScalar f0[]) {
        const PylithInt _dim = 3;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit_source_context(
            _dim, &poroelasticContext, &rheologyContext, f0);

    } // f0p_implicit_source

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static
    void f0p_implicit_source_body(const PylithInt dim,
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
                                  const PylithScalar x[],
                                  const PylithInt numConstants,
                                  const PylithScalar constants[],
                                  PylithScalar f0[]) {
        const PylithInt _dim = 3;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextBodyForceSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit_source_body_context(
            _dim, &poroelasticContext, &rheologyContext, f0);
    } // f0p_implicit_source_body

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_grav(const PylithInt dim,
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
                                  const PylithScalar x[],
                                  const PylithInt numConstants,
                                  const PylithScalar constants[],
                                  PylithScalar f0[]) {
        const PylithInt _dim = 3;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextGravitySourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit_source_grav_context(
            _dim, &poroelasticContext, &rheologyContext, f0);
       
    } // f0p_implicit_source_grav

        // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_grav_body(const PylithInt dim,
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
                                       const PylithScalar x[],
                                       const PylithInt numConstants,
                                       const PylithScalar constants[],
                                       PylithScalar f0[]) {
        const PylithInt _dim = 3;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextGravityBodyForceSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit_source_grav_body_context(
            _dim, &poroelasticContext, &rheologyContext, f0);
    } // f0p_implicit_source_grav_body

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear visco-poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Quasi - Static Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static inline
    void f1u(const PylithInt dim,
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
             const PylithScalar x[],
             const PylithInt numConstants,
             const PylithScalar constants[],
             PylithScalar f1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContextIsotropicLinearPoroelasticity;
        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContextIsotropicLinearMaxwell;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextQuasistatic(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContext(
            &rheologyContextIsotropicLinearMaxwell, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        Tensor strain;
        pylith::fekernels::Elasticity3D::infinitesimalStrain(strainContext, &strain);

        Tensor stress;
        TensorOps tensorOps = pylith::fekernels::Tensor::ops3D;
        pylith::fekernels::PoroIsotropicLinearMaxwell::cauchyStress(&rheologyContextIsotropicLinearPoroelasticity, &rheologyContextIsotropicLinearMaxwell, strain, tensorOps, &stress);

        PylithScalar stressTensor[9] = {0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
        tensorOps.toTensor(stress, stressTensor);

        for (PylithInt i = 0; i < _dim*_dim; ++i) {
                f1[i] -= stressTensor[i];
        } // for

    } // f1u

     // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static inline
    void f1u_refstate(const PylithInt dim,
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
                      const PylithScalar x[],
                      const PylithInt numConstants,
                      const PylithScalar constants[],
                      PylithScalar f1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContextIsotropicLinearPoroelasticity;
        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContextIsotropicLinearMaxwell;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextQuasistatic(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextRefState(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContext(
            &rheologyContextIsotropicLinearMaxwell, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContextRefState(
            &rheologyContextIsotropicLinearMaxwell, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        Tensor strain;
        pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain(strainContext, &strain);

        Tensor stress;
        TensorOps tensorOps = pylith::fekernels::Tensor::ops3D;
        pylith::fekernels::PoroIsotropicLinearMaxwell::cauchyStress_refState(&rheologyContextIsotropicLinearPoroelasticity, &rheologyContextIsotropicLinearMaxwell, strain, tensorOps, &stress);

        PylithScalar stressTensor[9] = {0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
        tensorOps.toTensor(stress, stressTensor);

        for (PylithInt i = 0; i < _dim*_dim; ++i) {
            f1[i] -= stressTensor[i];
        } // for


    } // f1u_refstate

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity
     *
     */
    static inline
    void f1p(const PylithInt dim,
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
             const PylithScalar x[],
             const PylithInt numConstants,
             const PylithScalar constants[],
             PylithScalar f1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_tensor_permeability(const PylithInt dim,
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
                                 const PylithScalar x[],
                                 const PylithInt numConstants,
                                 const PylithScalar constants[],
                                 PylithScalar f1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body(const PylithInt dim,
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
                  const PylithScalar x[],
                  const PylithInt numConstants,
                  const PylithScalar constants[],
                  PylithScalar f1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_body

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_tensor_permeability(const PylithInt dim,
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
                                      const PylithScalar x[],
                                      const PylithInt numConstants,
                                      const PylithScalar constants[],
                                      PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_body_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_gravity(const PylithInt dim,
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
                     const PylithScalar x[],
                     const PylithInt numConstants,
                     const PylithScalar constants[],
                     PylithScalar f1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_gravity

         // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_gravity_tensor_permeability(const PylithInt dim,
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
                                         const PylithScalar x[],
                                         const PylithInt numConstants,
                                         const PylithScalar constants[],
                                         PylithScalar f1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_gravity_tensor_permeability

     // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_gravity(const PylithInt dim,
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
                          const PylithScalar x[],
                          const PylithInt numConstants,
                          const PylithScalar constants[],
                          PylithScalar f1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravityBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_body_gravity

        // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_gravity_tensor_permeability(const PylithInt dim,
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
                                              const PylithScalar x[],
                                              const PylithInt numConstants,
                                              const PylithScalar constants[],
                                              PylithScalar f1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravityBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_body_gravity_tensor_permeability

        // =========================== LHS Jacobian ============================

    // ----------------------------------------------------------------------
    /* Jf3_uu entry function for isotropic linear visco-poroelasticity.
     */
    static inline
    void Jf3uu(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf3[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context isotropicLinearPoroelasticityRheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &isotropicLinearPoroelasticityRheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearMaxwell::Context isotropicLinearMaxwellRheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContext(
            &isotropicLinearMaxwellRheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        const PylithScalar shearModulus = isotropicLinearPoroelasticityRheologyContext.shearModulus;
        const PylithScalar maxwellTime = isotropicLinearMaxwellRheologyContext.maxwellTime;
        const PylithScalar dt = isotropicLinearMaxwellRheologyContext.dt;

        const PylithScalar dq = pylith::fekernels::IsotropicLinearMaxwell::viscousStrainCoeff(dt, maxwellTime);

        //Unique components of Jacobian.
        const PylithReal C1111 = 4.0/3.0 * shearModulus * dq;
        const PylithReal C1122 = -2.0/3.0 * shearModulus * dq;
        const PylithReal C1212 = shearModulus * dq;

                /* j(f,g,df,dg) = C(f,df,g,dg)
         *
         * 0:  j0000 = C1111 = bulkModulus + 4*dq*shearModulus/3
         * 1:  j0001 = C1112 = 0
         * 2:  j0002 = C1113 = 0
         * 3:  j0010 = C1211 = 0
         * 4:  j0011 = C1212 = dq*shearModulus
         * 5:  j0012 = C1213 = 0
         * 6:  j0020 = C1311 = 0
         * 7:  j0021 = C1312 = 0
         * 8:  j0022 = C1313 = dq*shearModulus
         * 9:  j0100 = C1121 = 0
         * 10:  j0101 = C1122 = bulkModulus - 2*dq*shearModulus/3
         * 11:  j0102 = C1123 = 0
         * 12:  j0110 = C1221 = dq*shearModulus
         * 13:  j0111 = C1222 = 0
         * 14:  j0112 = C1223 = 0
         * 15:  j0120 = C1321 = 0
         * 16:  j0121 = C1322 = 0
         * 17:  j0122 = C1323 = 0
         * 18:  j0200 = C1131 = 0
         * 19:  j0201 = C1132 = 0
         * 20:  j0202 = C1133 = bulkModulus - 2*dq*shearModulus/3
         * 21:  j0210 = C1231 = 0
         * 22:  j0211 = C1232 = 0
         * 23:  j0212 = C1233 = 0
         * 24:  j0220 = C1331 = dq*shearModulus
         * 25:  j0221 = C1332 = 0
         * 26:  j0222 = C1333 = 0
         * 27:  j1000 = C2111 = 0
         * 28:  j1001 = C2112 = dq*shearModulus
         * 29:  j1002 = C2113 = 0
         * 30:  j1010 = C2211 = bulkModulus - 2*dq*shearModulus/3
         * 31:  j1011 = C2212 = 0
         * 32:  j1012 = C2213 = 0
         * 33:  j1020 = C2311 = 0
         * 34:  j1021 = C2312 = 0
         * 35:  j1022 = C2313 = 0
         * 36:  j1100 = C2121 = dq*shearModulus
         * 37:  j1101 = C2122 = 0
         * 38:  j1102 = C2123 = 0
         * 39:  j1110 = C2221 = 0
         * 40:  j1111 = C2222 = bulkModulus + 4*dq*shearModulus/3
         * 41:  j1112 = C2223 = 0
         * 42:  j1120 = C2321 = 0
         * 43:  j1121 = C2322 = 0
         * 44:  j1122 = C2323 = dq*shearModulus
         * 45:  j1200 = C2131 = 0
         * 46:  j1201 = C2132 = 0
         * 47:  j1202 = C2133 = 0
         * 48:  j1210 = C2231 = 0
         * 49:  j1211 = C2232 = 0
         * 50:  j1212 = C2233 = bulkModulus - 2*dq*shearModulus/3
         * 51:  j1220 = C2331 = 0
         * 52:  j1221 = C2332 = dq*shearModulus
         * 53:  j1222 = C2333 = 0
         * 54:  j2000 = C3111 = 0
         * 55:  j2001 = C3112 = 0
         * 56:  j2002 = C3113 = dq*shearModulus
         * 57:  j2010 = C3211 = 0
         * 58:  j2011 = C3212 = 0
         * 59:  j2012 = C3213 = 0
         * 60:  j2020 = C3311 = bulkModulus - 2*dq*shearModulus/3
         * 61:  j2021 = C3312 = 0
         * 62:  j2022 = C3313 = 0
         * 63:  j2100 = C3121 = 0
         * 64:  j2101 = C3122 = 0
         * 65:  j2102 = C3123 = 0
         * 66:  j2110 = C3221 = 0
         * 67:  j2111 = C3222 = 0
         * 68:  j2112 = C3223 = dq*shearModulus
         * 69:  j2120 = C3321 = 0
         * 70:  j2121 = C3322 = bulkModulus - 2*dq*shearModulus/3
         * 71:  j2122 = C3323 = 0
         * 72:  j2200 = C3131 = dq*shearModulus
         * 73:  j2201 = C3132 = 0
         * 74:  j2202 = C3133 = 0
         * 75:  j2210 = C3231 = 0
         * 76:  j2211 = C3232 = dq*shearModulus
         * 77:  j2212 = C3233 = 0
         * 78:  j2220 = C3331 = 0
         * 79:  j2221 = C3332 = 0
         * 80:  j2222 = C3333 = bulkModulus + 4*dq*shearModulus/3
         */

        /* Nonzero Jacobian entries. */
        Jf3[0] -= C1111; /* j0000 */
        Jf3[4] -= C1212; /* j0011 */
        Jf3[8] -= C1212; /* j0022 */
        Jf3[10] -= C1122; /* j0101 */
        Jf3[12] -= C1212; /* j0110 */
        Jf3[20] -= C1122; /* j0202 */
        Jf3[24] -= C1212; /* j0220 */
        Jf3[28] -= C1212; /* j1001 */
        Jf3[30] -= C1122; /* j1010 */
        Jf3[36] -= C1212; /* j1100 */
        Jf3[40] -= C1111; /* j1111 */
        Jf3[44] -= C1212; /* j1122 */
        Jf3[50] -= C1122; /* j1212 */
        Jf3[52] -= C1212; /* j1221 */
        Jf3[56] -= C1212; /* j2002 */
        Jf3[60] -= C1122; /* j2020 */
        Jf3[68] -= C1212; /* j2112 */
        Jf3[70] -= C1122; /* j2121 */
        Jf3[72] -= C1212; /* j2200 */
        Jf3[76] -= C1212; /* j2211 */
        Jf3[80] -= C1111; /* j2222 */

    }

    // ----------------------------------------------------------------------
    /** Jf2_up entry function for isotropic linear poroviscoelasticity.
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf2up(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf2[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf2up_context(
            dim, &rheologyContext, Jf2);
    } // Jf2up

        // -----------------------------------------------------------------------------
    // Jf2ue function for isotropic linear poroviscoelasticity.
    static inline
    void Jf2ue(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf2[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf2ue_context(
            dim, &rheologyContext, Jf2);
    } // Jf2ue

        // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroviscoelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf3pp(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf3[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf3pp_context(
            dim, &rheologyContext, Jf3);

    } // Jf3pp

       // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroviscoelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf3pp_tensor_permeability(const PylithInt dim,
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
                                   const PylithReal s_tshift,
                                   const PylithScalar x[],
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar Jf3[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf3pp_context(
            dim, &rheologyContext, Jf3);

    } // Jf3pp_tensor_permeability

     // ----------------------------------------------------------------------
    /** Jf0_pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0pp(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf0[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf0pp_context(
            _dim, s_tshift, &rheologyContext, Jf0);
        
    } // Jf0pp

     // ----------------------------------------------------------------------
    /** Jf0_pe entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0pe(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf0[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf0pe_context(
            _dim, s_tshift, &rheologyContext, Jf0);

    } // Jf0pe

    // ===========================================================================================
    // Kernels for updating state variables
    // ===========================================================================================

    // Use pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector() to update total strain.

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating viscous strain as a vector for 2D plane strain isotropic
     * linear Maxwell viscoelasticity.
     *
     * Used to output viscous strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void viscousStrain_infinitesimalStrain_asVector(const PylithInt dim,
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
                                                    const PylithScalar x[],
                                                    const PylithInt numConstants,
                                                    const PylithScalar constants[],
                                                    PylithScalar viscousStrain[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearMaxwell::viscousStrain_asVector(
            strainContext, rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::Tensor::ops3D,
            viscousStrain);
    }

// ---------------------------------------------------------------------------------------------------------------------
    /* Update porosity for a linear poroelastic material, implicit.
     */
    static inline
    void updatePorosityImplicit(const PylithInt dim,
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
                                const PylithScalar x[],
                                const PylithInt numConstants,
                                const PylithScalar constants[],
                                PylithScalar porosity[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Incoming solution fields.
        const PylithInt i_pressure_t = 4;
        const PylithInt i_trace_strain_t = 5;

        // Incoming re-packed auxiliary field.

        // Poroelasticity
        const PylithInt i_porosityPrev = 3;

        // IsotropicLinearPoroelasticity
        const PylithInt i_drainedBulkModulus = numA - 7;
        const PylithInt i_biotCoefficient = numA - 6;

        // Constants
        const PylithScalar dt = constants[0];

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 6);
        assert(numA >= 12);
        assert(aOff);
        assert(aOff[i_porosityPrev] >= 0);
        assert(porosity);

        // Do stuff
        const PylithScalar pressure_t = s ? s[sOff[i_pressure_t]] : 0.0;
        const PylithScalar trace_strain_t = s ? s[sOff[i_trace_strain_t]] : 0.0;

        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar porosityPrev = a[aOff[i_porosityPrev]];

        // Update porosity
        porosity[0] = porosityPrev + dt * ((biotCoefficient - porosityPrev) * trace_strain_t +
                                           ((1.0 - biotCoefficient) * (biotCoefficient - porosityPrev)) /
                                           drainedBulkModulus * pressure_t);
        porosity[0] = std::max(0.0, std::min(1.0, porosity[0]));


    } // updatePorosityImplicit

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 2D isotropic linear poroelasticity with
     * infinitesimal strain WITHOUT reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., biot_coefficient(1), shear_modulus(1), drained_bulk_modulus(1)]
     */
    static inline
    void cauchyStress_infinitesimalStrain_asVector(const PylithInt dim,
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
                                                   const PylithScalar x[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar stressVector[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContextIsotropicLinearPoroelasticity;
        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContextIsotropicLinearMaxwell;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextQuasistatic(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContext(
            &rheologyContextIsotropicLinearMaxwell, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);


        Tensor strain;
        pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain(strainContext, &strain);

        Tensor stress;
        TensorOps tensorOps = pylith::fekernels::Tensor::ops3D;
        pylith::fekernels::PoroIsotropicLinearMaxwell::cauchyStress(&rheologyContextIsotropicLinearPoroelasticity, &rheologyContextIsotropicLinearMaxwell, strain, tensorOps, &stress);

        tensorOps.toVector(stress, stressVector);

    } // cauchyStress_infinitesimalStrain_asVector

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for DD isotropic linear poroelasticity with
     * infinitesimal strain WITH reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., biot_coefficient(1), shear_modulus(1), drained_bulk_modulus(1)]
     */
    static inline
    void cauchyStress_infinitesimalStrain_refState_asVector(const PylithInt dim,
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
                                                   const PylithScalar x[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar stressVector[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContextIsotropicLinearPoroelasticity;
        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContextIsotropicLinearMaxwell;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextRefState(
            &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        // Using dynamic formulation for trace strain, assuming that it will be equal to the variable
        // for QS
        // pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContextDynamic(
        //     &rheologyContextIsotropicLinearPoroelasticity, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
        //     t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContext(
            &rheologyContextIsotropicLinearMaxwell, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearMaxwellContextRefState(
            &rheologyContextIsotropicLinearMaxwell, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        Tensor strain;
        pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain(strainContext, &strain);

        Tensor stress;
        TensorOps tensorOps = pylith::fekernels::Tensor::ops3D;
        pylith::fekernels::PoroIsotropicLinearMaxwell::cauchyStress_refState(&rheologyContextIsotropicLinearPoroelasticity, &rheologyContextIsotropicLinearMaxwell, strain, tensorOps, &stress);

        tensorOps.toVector(stress, stressVector);

    } // cauchyStress_infinitesimalStrain_refState_asVector

        // Calculate water content
    static inline
    void waterContent_asScalar(const PylithInt dim,
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
                               const PylithScalar x[],
                               const PylithInt numConstants,
                               const PylithScalar constants[],
                               PylithReal* waterContent) {
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::PoroIsotropicLinearMaxwell::setIsotropicLinearPoroelasticityContext(
            &rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::waterContent_asScalar_context(
            dim, &poroelasticContext, &rheologyContext, waterContent);

    } // waterContent_asScalar

};// PoroIsotropicLinearMaxwell3D