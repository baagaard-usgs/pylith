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

#include "pylith/bc/bcfwd.hh" // forward declarations
#include "pylith/topology/SubfieldFactory.hh" // ISA SubfieldFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES TimeHistory

class pylith::bc::SubfieldFactory : public pylith::topology::SubfieldFactory {
    friend class TestDirichletAuxiliaryFactory; // unit testing

    // PUBLIC ENUMS ////////////////////////////////////////////////////////////////////////////////////////////////////
public:

    enum ComponentsType {
        XYZ=0, ///< Directions along coordinate axes.
        TANGENTIAL_NORMAL=1, ///< Directions tangential and normal to the boundary.
    };

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // Auxiliary subfields
    static const std::string density;
    static const std::string vp;
    static const std::string vs;
    static const std::string initial_amplitude;
    static const std::string rate_amplitude;
    static const std::string rate_start_time;
    static const std::string time_history_amplitude;
    static const std::string time_history_start_time;
    static const std::string time_history_value;

    // Diagnostic subfields
    static const std::string normal_dir;
    static const std::string tangential_dir_horiz;
    static const std::string tangential_dir_vert;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param[in] componentsType Reference for coordinate directions in auxiliary subfield.s
     */
    SubfieldFactory(const ComponentsType componentsType=XYZ);

    /// Destructor.
    ~SubfieldFactory(void) override;

    /** Set description of reference field.
     *
     * @param[in] refDescription Description of reference field.
     */
    void setRefDescription(const std::shared_ptr<pylith::topology::FieldBase::Description>& refDescription=nullptr);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Get subfield description.
     *
     * @param[in] subfieldName Name of subfield.
     * @returns Description of subfield.
     */
    pylith::topology::FieldBase::Description _getDescription(const std::string& subfieldName,
                                                             const size_t spaceDim) const override;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ComponentsType _componentsType; ///< Coordinate system reference for field components.

    /// Description of reference field.
    std::shared_ptr<pylith::topology::FieldBase::Description> _refDescription;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    SubfieldFactory(const SubfieldFactory &) = delete;
    const SubfieldFactory& operator=(const SubfieldFactory&) = delete;

}; // class SubfieldFactory

// End of file
