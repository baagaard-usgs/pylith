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

/*
 * The PETSc VTK viewer collects all fields and then writes them all
 * at once. This means we need to cache fields locally in order to
 * allow the output manager to reuse fields for dimensionalizing,
 * etc. Other writers do not suffer from this restriction, so we
 * implement this functionality in DataWriterVTK.
 */

#include "pylith/meshio/DataWriter.hh" // ISA DataWriter

#include "pylith/topology/topologyfwd.hh" // HOLDSA Field
#include "pylith/utils/petscfwd.h" // HASA PetscDM

// DataWriterVTK --------------------------------------------------------
/// Object for writing finite-element data to VTK file.
class pylith::meshio::DataWriterVTK : public DataWriter {
    friend class TestDataWriterVTKMesh; // unit testing
    friend class TestDataWriterVTKMaterial; // unit testing
    friend class TestDataWriterVTKSubmesh; // unit testing
    friend class TestDataWriterVTKBCMesh; // unit testing
    friend class TestDataWriterVTKFaultMesh; // unit testing
    friend class TestDataWriterVTKPoints; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    DataWriterVTK(void);

    /// Destructor
    ~DataWriterVTK(void) override;

    /// Deallocate PETSc and local data structures.
    void deallocate(void) override;

    /** Set filename for VTK file.
     *
     * @param[in] filename Name of VTK file.
     */
    void setFilename(const char* filename);

    /** Set time format for time stamp in name of VTK file.
     *
     * @param[in] format C style time format for filename.
     */
    void setTimeFormat(const char* format);

    /** Set value used to normalize time stamp in name of VTK file.
     *
     * Time stamp is divided by this value (time in seconds).
     *
     * @param[in] value Value (time in seconds) used to normalize time stamp in
     * filename.
     */
    void setTimeConstant(const pylith::real value);

    /** Set precision of floating point values in output.
     *
     * @param[in] value Precision for floating point values.
     */
    void setPrecision(const int value);

    /** Prepare for writing files.
     *
     * @param[in] mesh Finite-element mesh.
     * @param[in] isInfo True if only writing info values.
     */
    void open(const topology::Mesh& mesh,
              const bool isInfo) override;

    /// Close output files.
    void close(void) override;

    /** Prepare file for data at a new time step.
     *
     * @param[in] t Time stamp for new data
     * @param[in] mesh Finite-element mesh.
     */
    void openTimeStep(const pylith::real t,
                      const topology::Mesh& mesh) override;

    /// Cleanup after writing data for a time step.
    void closeTimeStep(void) override;

    /** Write field over vertices to file.
     *
     * @param[in] t Time associated with field.
     * @param[in] subfield Subfield with basis order 1.
     */
    void writeVertexField(const pylith::real t,
                          const pylith::meshio::OutputSubfield& field) override;

    /** Write field over cells to file.
     *
     * @param[in] t Time associated with field.
     * @param[in] subfield Subfield with basis order 0.
     */
    void writeCellField(const pylith::real t,
                        const pylith::meshio::OutputSubfield& subfield) override;

    // PRIVATE METHODS //////////////////////////////////////////////////////
private:

    /** Generate filename for VTK file.
     *
     * @param[in] t Time in seconds.
     */
    std::string _getVTKFilename(const pylith::real t) const;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    DataWriterVTK(const DataWriterVTK&) = delete;
    const DataWriterVTK& operator=(const DataWriterVTK&) = delete;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    /// Time value (in seconds) used to normalize time stamp.
    pylith::real _timeConstant;
    int _precision; ///< Precision of floating point values in output.

    std::string _filename; ///< Name of VTK file.
    std::string _timeFormat; ///< C style time format for time stamp.

    PetscViewer _viewer; ///< Output file
    PetscDM _dm; ///< Handle to PETSc DM for mesh

    bool _isOpenTimeStep; ///< true if called openTimeStep().
    bool _wroteVertexHeader; ///< True if wrote header for vertex data.
    bool _wroteCellHeader; ///< True if wrote header for cell data

}; // DataWriterVTK

// End of file
