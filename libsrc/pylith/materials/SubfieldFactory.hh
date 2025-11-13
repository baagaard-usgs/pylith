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

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // Auxiliary subfields
    static const std::string density;
    static const std::string body_force;
    static const std::string gravity_field;
    static const std::string shear_modulus;
    static const std::string bulk_modulus;
    static const std::string reference_stress;
    static const std::string reference_strain;

    static const std::string maxwell_time;
    static const std::string genmaxwell_maxwell_time;
    static const std::string genmaxwell_shear_modulus_ratio;
    static const std::string powerlaw_reference_strain_rate;
    static const std::string powerlaw_reference_stress;
    static const std::string powerlaw_exponent;
    static const std::string total_strain;
    static const std::string deviatoric_stress;
    static const std::string viscous_strain;
    static const std::string genmaxwell_viscous_strain;

    static const std::string isotropic_permeability;
    static const std::string tensor_permeability;
    static const std::string drained_bulk_modulus;
    static const std::string biot_coefficient;
    static const std::string biot_modulus;
    static const std::string porosity;
    static const std::string solid_density;
    static const std::string fluid_density;
    static const std::string fluid_viscosity;
    static const std::string source_density;

    // Derived subfields
    static const std::string cauchy_stress;
    static const std::string cauchy_strain;
    static const std::string bulk_density;
    static const std::string water_content;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    SubfieldFactory(void);

    /// Destructor.
    ~SubfieldFactory(void) override;

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Get subfield description.
     *
     * @param[in] subfieldName Name of subfield.
     * @returns Description of subfield.
     */
    pylith::topology::FieldBase::Description _getDescription(const std::string& subfieldName,
                                                             const size_t spaceDim) const override;

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    SubfieldFactory(const SubfieldFactory &) = delete;
    const SubfieldFactory& operator=(const SubfieldFactory&) = delete;

}; // class SubfieldFactory

// End of file
