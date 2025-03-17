// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/materials/materialsfwd.hh" // forward declarations
#include "pylith/topology/SubfieldFactory.hh" // ISA SubfieldFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

class pylith::materials::SubfieldFactory : public pylith::topology::SubfieldFactory {
    friend class TestSubfieldFactory; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    SubfieldFactory(void);

    /// Destructor.
    virtual ~SubfieldFactory(void);

    /// Add density subfield to auxiliary subfields.
    void addDensity(void);

    /// Add body force subfield to auxiliary subfields.
    void addBodyForce(void);

    /** Add gravity subfield to auxiliary subfields.
     *
     * @param[in] gravityField Gravity field.
     */
    void addGravityField(std::shared_ptr<spatialdata::spatialdb::GravityField>& gravityField);

    /// Add shear modulus subfield to auxiliary subfields.
    void addShearModulus(void);

    /// Add bulk subfield to auxiliary subfields.
    void addBulkModulus(void);

    /// Add reference stress subfield to auxiliary subfields.
    void addReferenceStress(void);

    /// Add reference strain subfield to auxiliary subfields.
    void addReferenceStrain(void);

    /// Add Cauchy stress subfield to derived field.
    void addCauchyStress(void);

    /// Add Cauchy (infinitesimal) strain subfield to derived field.
    void addCauchyStrain(void);

    /// Add Maxwell time subfield to auxiliary subfields.
    void addMaxwellTime(void);

    /// Add Maxwell time subfield for Generalized Maxwell to auxiliary subfields.
    void addMaxwellTimeGeneralizedMaxwell(void);

    /// Add shear modulus ratio subfield for Generalized Maxwell to auxiliary subfields.
    void addShearModulusRatioGeneralizedMaxwell(void);

    /// Add power-law reference strain rate subfield to auxiliary subfields.
    void addPowerLawReferenceStrainRate(void);

    /// Add power-law reference stress subfield to auxiliary subfields.
    void addPowerLawReferenceStress(void);

    /// Add power-law exponent subfield to auxiliary subfields.
    void addPowerLawExponent(void);

    /// Add total strain subfield to auxiliary subfields.
    void addTotalStrain(void);

    /// Add stress subfield to auxiliary subfields.
    void addDeviatoricStress(void);

    /// Add viscous strain subfield to auxiliary subfields.
    void addViscousStrain(void);

    /// Add viscous strain subfield for Generalized Maxwell to auxiliary subfields.
    void addViscousStrainGeneralizedMaxwell(void);

    /// Add porosity subfield to auxiliary subfields.
    void addPorosity(void);

    /// Add solid density subfield to auxiliary subfields.
    void addSolidDensity(void);

    /// Add fluid density subfield to auxiliary subfields.
    void addFluidDensity(void);

    /// Add fluid viscosity subfield to auxiliary subfields.
    void addFluidViscosity(void);

    /// Add reference sourceDensity subfield to auxiliary fields.
    void addSourceDensity(void);

    /// Add isotropic permeability subfield to auxiliary subfields.
    void addIsotropicPermeability(void);

    /// Add tensor permeability subfield to auxiliary subfields.
    void addTensorPermeability(void);

    /// Add drained Bulk Modulus subfield to auxiliary subfields.
    void addDrainedBulkModulus(void);

    /// Add fluid Biot Coefficient subfield to auxiliary subfields.
    void addBiotCoefficient(void);

    /// Add fluid Biot Modulus subfield to auxiliary subfields.
    void addBiotModulus(void);

    /// Add bulk density subfield to derived field
    void addBulkDensity(void);

    /// Add water content subfield to derived field
    void addWaterContent(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    SubfieldFactory(const SubfieldFactory &); ///< Not implemented.
    const SubfieldFactory& operator=(const SubfieldFactory&); ///< Not implemented

}; // class SubfieldFactory

// End of file
