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
 * HDF5 schema for PyLith output.
 *
 * / - root group
 *   geometry - group
 *     coordsys - attribute string with coordinate system
 *     vertices - dataset [nvertices, spacedim]
 *   topology - group
 *     cell_type - attribute string with cell type
 *     cells - dataset [ncells, ncorners]
 *   vertex_fields - group
 *     VERTEX_FIELD (name of vertex field) - dataset
 *       [ntimesteps, nvertices, fiberdim]
 *   cell_fields - group
 *     CELL_FIELD (name of cell field) - dataset
 *       [ntimesteps, ncells, fiberdim]
 */

#include "pylith/meshio/DataWriter.hh" // ISA DataWriter

#include <string> // USES std::string
#include <map> // HASA std::map

class pylith::meshio::DataWriterHDF5Ext : public DataWriter { // DataWriterHDF5Ext
    friend class TestDataWriterHDF5ExtMesh; // unit testing
    friend class TestDataWriterHDF5ExtSubmesh; // unit testing
    friend class TestDataWriterHDF5ExtPoints; // unit testing
    friend class TestDataWriterHDF5ExtBCMesh; // unit testing
    friend class TestDataWriterHDF5ExtFaultMesh; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    DataWriterHDF5Ext(void);

    /// Destructor
    ~DataWriterHDF5Ext(void) override;

    /// Deallocate PETSc and local data structures.
    void deallocate(void) override;

    /** Set filename for HDF5 file.
     *
     * @param[in] filename Name of HDF5 file.
     */
    void setFilename(const char* filename);

    /** Generate filename for HDF5 file.
     *
     * Appends _info if only writing parameters.
     *
     * :KLUDGE: We should separate generating "info" files from the
     * DataWriter interface.
     *
     * @returns String for HDF5 filename.
     */
    std::string getHDF5Filename(void) const;

    /** Prepare for writing files.
     *
     * @param[in] mesh Finite-element mesh.
     * @param[in] isInfo True if only writing info values. */
    void open(const topology::Mesh& mesh,
              const bool isInfo) override;

    /// Close output files.
    void close(void) override;

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

    /** Write dataset with names of points to file.
     *
     * @param[in] names Array with name for each point, e.g., station name.
     * @param[in] mesh Finite-element mesh.
     *
     * Primarily used with OutputSolnPoints.
     */
    void writePointNames(const pylith::string_vector& names,
                         const topology::Mesh& mesh) override;

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    /// Generate filename for external dataset file.
    std::string _getDatasetFilename(const char* field) const;

    /** Write time stamp to file.
     *
     * @param[in] t Time in seconds.
     */
    void _writeTimeStamp(const pylith::real t);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    DataWriterHDF5Ext(const DataWriterHDF5Ext&) = delete;
    const DataWriterHDF5Ext& operator=(const DataWriterHDF5Ext&) = delete;

    // PRIVATE STRUCTS ////////////////////////////////////////////////////////////////////////////
private:

    struct ExternalDataset {
        PetscViewer viewer;
        pylith::integer numTimeSteps;
        pylith::integer numPoints;
        pylith::integer fiberDim;
    };
    typedef std::map<std::string, ExternalDataset> dataset_t;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    std::string _filename; ///< Name of HDF5 file.
    std::unique_ptr<HDF5> _h5; ///< HDF5 file
    dataset_t _datasets; ///< Datasets
    size_t _tstampIndex; ///< Index of last time stamp written.

}; // DataWriterHDF5Ext

// End of file
