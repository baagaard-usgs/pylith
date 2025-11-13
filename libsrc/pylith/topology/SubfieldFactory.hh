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

#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Discretization

#include "spatialdata/units/unitsfwd.hh" // HOLDSA Normalizer

#include <memory> // HASA std::shared_ptr

class pylith::topology::SubfieldFactory : public pylith::utils::GenericComponent {
    friend class TestSubfieldFactory; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    SubfieldFactory(void);

    /// Destructor.
    virtual ~SubfieldFactory(void);

    /** Set discretization information for subfield.
     *
     * @param[in] subfieldName Name of subfield.
     * @param[in] discretization Discretization for subfield.
     */
    void setDiscretization(const char* subfieldName,
                           const pylith::topology::FieldBase::Discretization& discretization);

    /** Check if factory has discretization information for subfield.
     *
     * @param[in] subfieldName Name of subfield.
     * @returns True if factory contains discretization for subfield.
     */
    bool hasDiscretization(const char* subfieldName);

    /** Open factory for setting up subfields.
     *
     * @param[inout] field Field for which subfields are to be created.
     * @param[in] normalizer Scales for nondimensionalization.
     * @param[in] spaceDim Spatial dimension of problem.
     */
    virtual
    void open(std::shared_ptr<pylith::topology::Field>& field,
              const std::shared_ptr<spatialdata::units::Nondimensional>& normalizer);

    /** Add subfield to field.
     *
     * @param[in] subfieldName Name of subfield.
     */
    void addSubfield(const std::string& subfieldName);

    /// Close factory after setting up subfields.
    void close(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Get subfield description.
     *
     * @param[in] subfieldName Name of subfield.
     * @param[in] spaceDim Spatial dimension.
     * @returns Description of subfield.
     */
    virtual
    pylith::topology::FieldBase::Description _getDescription(const std::string& subfieldName,
                                                             const size_t spaceDim) const = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    typedef std::map<std::string, pylith::topology::FieldBase::Discretization> subfields_t;

    std::shared_ptr<pylith::topology::Field> _field; ///< Field.
    subfields_t _subfields; ///< Discretization for each subfield.

    std::shared_ptr<spatialdata::units::Nondimensional> _normalizer; ///< Scales for nondimensionalization.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    SubfieldFactory(const SubfieldFactory &) = delete;
    const SubfieldFactory& operator=(const SubfieldFactory&) = delete;

}; // class SubfieldFactory

// End of file
