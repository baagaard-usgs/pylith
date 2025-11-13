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

#include "pylith/meshio/meshiofwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/array.hh" // HASA string_vector
#include <map> // HASA std::map

class pylith::meshio::OutputObserver : public pylith::utils::PyreComponent {
    friend class TestOutputObserver; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    OutputObserver(void);

    /// Destructor
    virtual ~OutputObserver(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set output trigger for how often to write output.
     *
     * @param[in] trigger Output trigger.
     */
    void setTrigger(const std::shared_ptr<pylith::meshio::OutputTrigger>& trigger);

    /** Get trigger for how often to write otuput.
     *
     * @returns Output trigger.
     */
    const std::shared_ptr<pylith::meshio::OutputTrigger>& getTrigger(void) const;

    /** Set writer to write data to file.
     *
     * @param[in] datawriter Writer for data.
     */
    void setWriter(const std::shared_ptr<pylith::meshio::DataWriter>& writer);

    /** Set basis order for output.
     *
     * @param[in] value Basis order for output.
     */
    void setOutputBasisOrder(const size_t value);

    /** Set number of mesh refinement levels for output.
     *
     * @param[in] value Number of mesh refinement levels for output.
     */
    void setRefineLevels(const size_t value);

    /** Set time scale.
     *
     * @param[in] value Time scale for dimensionalizing time.
     */
    virtual
    void setTimeScale(const pylith::real value);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Set context.
     *
     * @param[in] mesh Mesh associated with output.
     */
    void _setContext(const pylith::topology::Mesh& mesh);

    /** Get mesh associated with subfield output.
     *
     * @param[in] subfield Subfield for output.
     * @returns Mesh associated with output.
     */
    pylith::topology::Mesh* const _getOutputMesh(const pylith::meshio::OutputSubfield& subfield);

    /** Get output subfield, creating if necessary.
     *
     * @param[in] field Field containing subfields.
     * @param[in] submesh Submesh associated with output.
     * @param[in] name Name of subfield.
     */
    virtual
    OutputSubfield* const _getSubfield(const pylith::topology::Field& field,
                                       const pylith::topology::Mesh& submesh,
                                       const char* name);

    /** Append subfield at current time to output.
     *
     * @param[in] t Current time.
     * @param[in] subfield Subfield to write.
     */
    void _appendField(const pylith::real t,
                      const pylith::meshio::OutputSubfield& subfield);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    typedef std::map<std::string, std::unique_ptr<OutputSubfield> > subfield_t;

    pylith::real _timeScale; ///< Time scale for dimentionalizing time.
    subfield_t _subfields; ///< Subfields extracted for output.
    std::unique_ptr<pylith::topology::Mesh> _outputMesh; ///< Mesh associated with output.ß
    std::shared_ptr<DataWriter> _writer; ///< Writer for data.
    std::shared_ptr<OutputTrigger> _trigger; ///< Trigger for deciding how often to write output.
    pylith::integer _outputBasisOrder; ///< Basis order for output.
    pylith::integer _refineLevels; ///< Number of mesh refinement levels for output.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    OutputObserver(const OutputObserver&) = delete;
    const OutputObserver& operator=(const OutputObserver&) = delete;

}; // OutputObserver

// End of file
