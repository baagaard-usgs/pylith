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

#include "TestDataWriterHDF5.hh" // Implementation of class methods

#include "pylith/meshio/HDF5.hh" // USES HDF5
#include "pylith/utils/types.hh" // HASA PylithScalar
#include "pylith/utils/error.hh" // HASA PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <hdf5.h> // USES HDF5 API

#include <memory> // USES std::unique_ptr
#include <cstdint>
#include <map>
#include <numeric>

#if H5_VERSION_GE(1,12,0)
#define PYLITH_HDF5_USE_API_112
#endif

namespace pylith {
    namespace meshio {
        namespace _TestDataWriterHDF5 {
            class MeshData {
public:

                std::string name;
                size_t cellDim; // :TODO: Set and check
                std::vector<PylithReal> vertices;
                std::vector<PylithInt> cells;
                std::vector<PylithReal> time;

                struct Field {
                    std::vector<PylithScalar> data;
                    size_t numTimes;
                    size_t numComponents;
                    size_t numPoints;
                    std::string vectorFieldType; // :TODO: Set and check
                };
                std::map<std::string, Field> vertexFields;
                std::map<std::string, Field> cellFields;

                std::vector<size_t> cellsOrder; // :TODO: Set and use

                size_t numVertices;
                size_t spaceDim;
                size_t numCells;
                size_t numCorners;

public:

                MeshData(void);

                /** Load data from HDF5 file. */
                static
                MeshData* load(const char* filename);

                /** Verify that candidate data matches reference data. */
                static
                void verify(const MeshData& reference,
                            const MeshData& candidate);


private:

                void _loadGeometry(pylith::meshio::HDF5& h5);

                void _loadTopology(pylith::meshio::HDF5& h5);

                void _loadTime(pylith::meshio::HDF5& h5);

                static
                std::map<std::string, Field> _loadFields(pylith::meshio::HDF5& h5,
                                                         const char* group);

                std::vector<size_t> _getVerticesOrder(void) const;

                std::vector<size_t> _getCellsOrder(void) const;

                std::vector<size_t> _getCellCornersOrder(const size_t cell) const;

                static void _checkGeometry(const MeshData& reference,
                                           const MeshData& candidate);

                static void _checkTopology(const MeshData& reference,
                                           const MeshData& candidate);

                static void _checkTime(const MeshData& reference,
                                       const MeshData& candidate);

                static void _checkVertexFields(const MeshData& reference,
                                               const MeshData& candidate);

                static void _checkCellFields(const MeshData& reference,
                                             const MeshData& candidate);


            }; // MeshData
        } // _TestDataWriterHDF5

    }
}


// ------------------------------------------------------------------------------------------------
pylith::meshio::_TestDataWriterHDF5::MeshData::MeshData(void) {}


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_TestDataWriterHDF5::MeshData::_loadGeometry(pylith::meshio::HDF5& h5) {
    PYLITH_METHOD_BEGIN;

    REQUIRE(h5.hasDataset("/geometry/vertices"));
    REQUIRE(sizeof(PylithReal) == sizeof(double));

    HDF5::Dataset<double> dataset = h5.readDataset<double>("/geometry", "vertices");
    REQUIRE(dataset.dims.size() == size_t(2));
    numVertices = dataset.dims[0];
    spaceDim = dataset.dims[1];
    vertices = dataset.data;

    PYLITH_METHOD_END;
} // _loadGeometry


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_TestDataWriterHDF5::MeshData::_loadTopology(pylith::meshio::HDF5& h5) {
    PYLITH_METHOD_BEGIN;

    REQUIRE(h5.hasDataset("/viz/topology/cells"));
    REQUIRE(sizeof(size_t) == sizeof(uint64_t));

    HDF5::Dataset<PylithInt> dataset = h5.readDataset<PylithInt>("/viz/topology", "cells");
    REQUIRE(dataset.dims.size() == 2);
    numCells = dataset.dims[0];
    numCorners = dataset.dims[1];
    cells = dataset.data;

    PYLITH_METHOD_END;
} // _loadTopology


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_TestDataWriterHDF5::MeshData::_loadTime(pylith::meshio::HDF5& h5) {
    PYLITH_METHOD_BEGIN;

    if (h5.hasDataset("/time")) {
        REQUIRE(sizeof(PylithReal) == sizeof(double));

        HDF5::Dataset<double> dataset = h5.readDataset<double>("/", "time");
        REQUIRE(dataset.dims.size() == 3);
        CHECK(dataset.dims[1] == 1);
        CHECK(dataset.dims[2] == 1);
        time = dataset.data;
    } else {
        time.clear();
    } // if/else

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
std::map<std::string, pylith::meshio::_TestDataWriterHDF5::MeshData::Field>
pylith::meshio::_TestDataWriterHDF5::MeshData::_loadFields(pylith::meshio::HDF5& h5,
                                                           const char* group) {
    PYLITH_METHOD_BEGIN;

    std::map<std::string, Field> fields;

    if (!h5.hasGroup(group)) {
        PYLITH_METHOD_RETURN(fields);
    } // if

    pylith::string_vector names = h5.getGroupDatasets(group);

    for (const std::string& name : names) {
        HDF5::Dataset<double> dataset = h5.readDataset<double>(group, name.c_str());

        Field field;
        field.numTimes = dataset.dims[0];
        field.numPoints = dataset.dims[1];
        field.numComponents = dataset.dims[2];

        field.data = dataset.data;

        std::string fieldPath = std::string(group) + std::string("/") + name;
        field.vectorFieldType = h5.readAttributeString(fieldPath.c_str(), "vector_field_type");

        fields.emplace(name, field);
    } // for

    PYLITH_METHOD_RETURN(fields);
}


// ------------------------------------------------------------------------------------------------
std::vector<size_t>
pylith::meshio::_TestDataWriterHDF5::MeshData::_getVerticesOrder(void) const {
    std::vector<size_t> indices(numVertices);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](size_t iVertex, size_t jVertex) {
        for (size_t iDim = 0; iDim < spaceDim; ++iDim) {
            double iCoord = vertices[iVertex * spaceDim + iDim];
            double jCoord = vertices[jVertex * spaceDim + iDim];
            if (iCoord != jCoord) {return iCoord < jCoord;}
        }
        return false;
    }); // sort
    return indices;
} // _getVerticesOder


// ------------------------------------------------------------------------------------------------
std::vector<size_t>
pylith::meshio::_TestDataWriterHDF5::MeshData::_getCellsOrder(void) const {
    // Compute cell centroids
    std::vector<PylithReal> centroids(numCells*spaceDim, 0.0);
    for (size_t iCell = 0; iCell < numCells; ++iCell) {
        for (size_t iCorner = 0; iCorner < numCorners; ++iCorner) {
            size_t indexVertex = cells[iCell * numCorners + iCorner];
            for (size_t iDim = 0; iDim < spaceDim; ++iDim) {
                centroids[iCell*spaceDim+iDim] += vertices[indexVertex * spaceDim + iDim];
            } // for
        } // for
    } // for
    for (double& x : centroids) {
        x /= PylithReal(numCorners);
    } // for

    std::vector<size_t> indices(numCells);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](size_t iCell, size_t jCell) {
        for (hsize_t iDim = 0; iDim < spaceDim; ++iDim) {
            double iCoord = centroids[iCell * spaceDim + iDim];
            double jCoord = centroids[jCell * spaceDim + iDim];
            if (iCoord != jCoord) {return iCoord < jCoord;}
        }
        return false;
    }); // sort
    return indices;
} // _getCellsOrder


// ------------------------------------------------------------------------------------------------
std::vector<size_t>
pylith::meshio::_TestDataWriterHDF5::MeshData::_getCellCornersOrder(const size_t cellIndex) const {
    REQUIRE(cellIndex < numCells);

    std::vector<PylithReal> cellVertices(numCorners*spaceDim);
    for (size_t iCorner = 0; iCorner < numCorners; ++iCorner) {
        const size_t vertexIndex = cells[cellIndex*numCorners+iCorner];
        for (size_t iDim = 0; iDim < spaceDim; ++iDim) {
            cellVertices[iCorner*spaceDim+iDim] = vertices[vertexIndex*spaceDim+iDim];
        } // for
    } // for

    std::vector<size_t> indicesV(numCorners);
    std::iota(indicesV.begin(), indicesV.end(), 0);
    std::sort(indicesV.begin(), indicesV.end(), [&](size_t iVertex, size_t jVertex) {
        for (size_t iDim = 0; iDim < spaceDim; ++iDim) {
            double iCoord = cellVertices[iVertex * spaceDim + iDim];
            double jCoord = cellVertices[jVertex * spaceDim + iDim];
            if (iCoord != jCoord) {return iCoord < jCoord;}
        }
        return false;
    }); // sort

    return indicesV;
} // _getCellCornersOrder


// ------------------------------------------------------------------------------------------------
// Load data from HDF5 file.
pylith::meshio::_TestDataWriterHDF5::MeshData*
pylith::meshio::_TestDataWriterHDF5::MeshData::load(const char* filename) {
    PYLITH_METHOD_BEGIN;

    MeshData* data = new MeshData;
    try {
        pylith::meshio::HDF5 h5(filename, H5F_ACC_RDONLY);
        data->name = filename;
        data->_loadGeometry(h5);
        data->_loadTopology(h5);
        data->_loadTime(h5);
        data->vertexFields = _loadFields(h5, "/vertex_fields");
        data->cellFields = _loadFields(h5, "/cell_fields");
        h5.close();
    } catch (const std::exception& err) {
        delete data;
        std::ostringstream msg;
        msg << "I/O error loading data from " << filename << "." << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        delete data;
        std::ostringstream msg;
        msg << "Unknown I/O error while loading data from " << filename << ".";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_RETURN(data);
} // load


// ------------------------------------------------------------------------------------------------
// Verify that candidate data matches reference data.
void
pylith::meshio::_TestDataWriterHDF5::MeshData::verify(const MeshData& reference,
                                                      const MeshData& candidate) {
    PYLITH_METHOD_BEGIN;

    INFO("Checking candidate '" << candidate.name << "' against reference '" << reference.name << "'.");

    _checkGeometry(reference, candidate);
    _checkTopology(reference, candidate);
    _checkTime(reference, candidate);
    _checkVertexFields(reference, candidate);
    _checkCellFields(reference, candidate);

    PYLITH_METHOD_END;
} // verify


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_TestDataWriterHDF5::MeshData::_checkGeometry(const MeshData& reference,
                                                              const MeshData& candidate) {
    PYLITH_METHOD_BEGIN;

    INFO("Checking geometry.");
    REQUIRE(reference.numVertices == candidate.numVertices);
    REQUIRE(reference.spaceDim == candidate.spaceDim);

    const size_t spaceDim = reference.spaceDim;
    const std::vector<size_t>& referenceOrder = reference._getVerticesOrder();
    const std::vector<size_t>& candidateOrder = candidate._getVerticesOrder();
    const PylithScalar tolerance = 1.0e-6;
    for (size_t iVertex = 0; iVertex < reference.numVertices; ++iVertex) {
        const size_t vertexReference = referenceOrder[iVertex];
        const size_t vertexCandidate = candidateOrder[iVertex];
        for (size_t iDim = 0; iDim < spaceDim; ++iDim) {
            const size_t indexC = vertexCandidate*spaceDim+iDim;
            const size_t indexR = vertexReference*spaceDim+iDim;
            const PylithScalar toleranceV = std::max(tolerance, tolerance*reference.vertices[indexR]);
            INFO("Checking coordinates of candidate vertex " << vertexCandidate << " against reference vertex " << vertexReference << " for dim " << iDim << ".");
            CHECK_THAT(candidate.vertices[indexC], Catch::Matchers::WithinAbs(reference.vertices[indexR], toleranceV));
        } // for
    } // for

    PYLITH_METHOD_END;
} // _checkGeometry


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_TestDataWriterHDF5::MeshData::_checkTopology(const MeshData& reference,
                                                              const MeshData& candidate) {
    PYLITH_METHOD_BEGIN;

    INFO("Checking topology.");
    REQUIRE(reference.numCells == candidate.numCells);
    REQUIRE(reference.numCorners == candidate.numCorners);

    // Geometry information used in checking topology
    REQUIRE(reference.numVertices == candidate.numVertices);
    REQUIRE(reference.spaceDim == candidate.spaceDim);

    const size_t spaceDim = reference.spaceDim;
    const size_t numCorners = reference.numCorners;
    const std::vector<size_t>& referenceOrder = reference._getCellsOrder();
    const std::vector<size_t>& candidateOrder = candidate._getCellsOrder();

    const PylithScalar tolerance = 1.0e-6;
    for (size_t iCell = 0; iCell < reference.numCells; ++iCell) {
        const size_t cellReference = referenceOrder[iCell];
        const size_t cellCandidate = candidateOrder[iCell];
        INFO("Checking candidate cell " << cellCandidate << " against reference cell " << cellReference << ".");
        const std::vector<size_t>& cellVerticesOrderR = reference._getCellCornersOrder(referenceOrder[iCell]);
        const std::vector<size_t>& cellVerticesOrderC = candidate._getCellCornersOrder(candidateOrder[iCell]);

        for (size_t iCorner = 0; iCorner < numCorners; ++iCorner) {
            const size_t iVertexC = candidate.cells[cellCandidate*numCorners+cellVerticesOrderC[iCorner]];
            const size_t iVertexR = reference.cells[cellReference*numCorners+cellVerticesOrderR[iCorner]];
            INFO("Checking candidate vertex " << iVertexC << " against reference vertex " << iVertexR << ".");
            for (size_t iDim = 0; iDim < spaceDim; ++iDim) {
                const PylithScalar toleranceV = std::max(tolerance, tolerance*reference.vertices[iVertexR]);
                CHECK_THAT(candidate.vertices[iVertexC*spaceDim+iDim], Catch::Matchers::WithinAbs(reference.vertices[iVertexR*spaceDim+iDim], toleranceV));
            } // for
        } // for
    } // for

    PYLITH_METHOD_END;
} // _checkTopology


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_TestDataWriterHDF5::MeshData::_checkTime(const MeshData& reference,
                                                          const MeshData& candidate) {
    PYLITH_METHOD_BEGIN;

    INFO("Checking time field.");

    REQUIRE(reference.time.size() == candidate.time.size());
    const size_t numTimes = reference.time.size();

    const PylithScalar tolerance = 1.0e-6;
    for (size_t iTime = 0; iTime < numTimes; ++iTime) {
        INFO("Checking time " << iTime << ".");
        const PylithScalar toleranceV = std::max(tolerance, tolerance*reference.time[iTime]);
        CHECK_THAT(candidate.time[iTime], Catch::Matchers::WithinAbs(reference.time[iTime], toleranceV));
    } // for

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_TestDataWriterHDF5::MeshData::_checkVertexFields(const MeshData& reference,
                                                                  const MeshData& candidate) {
    PYLITH_METHOD_BEGIN;

#if 0
    REQUIRE(candidate.vertexFields.size() == reference.vertexFields.size());
#endif
    for (const auto& kv : reference.vertexFields) {
        const std::string& name = kv.first;
        const Field& referenceField = kv.second;

        const auto& iterC = candidate.vertexFields.find(name);
        REQUIRE(iterC != candidate.vertexFields.end());
        const Field& candidateField = iterC->second;

        INFO("Checking vertex field '" << name << "'.");

        REQUIRE(reference.numVertices == candidate.numVertices);
        REQUIRE(reference.spaceDim == candidate.spaceDim);
        REQUIRE(reference.numCells == candidate.numCells);
        REQUIRE(reference.numCorners == candidate.numCorners);
        REQUIRE(referenceField.numPoints == candidateField.numPoints);
        REQUIRE(referenceField.numComponents == candidateField.numComponents);

        const size_t numComponents = referenceField.numComponents;
        const size_t spaceDim = reference.spaceDim;
        const size_t numCorners = reference.numCorners;
        const std::vector<size_t>& referenceOrder = reference._getCellsOrder();
        const std::vector<size_t>& candidateOrder = candidate._getCellsOrder();

        const PylithScalar tolerance = 1.0e-6;
        for (size_t iCell = 0; iCell < reference.numCells; ++iCell) {
            const size_t cellReference = referenceOrder[iCell];
            const size_t cellCandidate = candidateOrder[iCell];
            INFO("Checking candidate cell " << cellCandidate << " against reference cell " << cellReference << ".");
            const std::vector<size_t>& cellVerticesOrderR = reference._getCellCornersOrder(referenceOrder[iCell]);
            const std::vector<size_t>& cellVerticesOrderC = candidate._getCellCornersOrder(candidateOrder[iCell]);

            for (size_t iCorner = 0; iCorner < numCorners; ++iCorner) {
                const size_t iVertexC = candidate.cells[cellCandidate*numCorners+cellVerticesOrderC[iCorner]];
                const size_t iVertexR = reference.cells[cellReference*numCorners+cellVerticesOrderR[iCorner]];
                INFO("Checking candidate vertex " << iVertexC << " against reference vertex " << iVertexR << ".");
                for (size_t iComponent = 0; iComponent < numComponents; ++iComponent) {
                    const PylithScalar toleranceV = std::max(tolerance, tolerance*referenceField.data[iVertexR*numComponents+iComponent]);
#if 0
                    CHECK_THAT(candidateField.data[iVertexC*numComponents+iComponent],
                               Catch::Matchers::WithinAbs(referenceField.data[iVertexR*numComponents+iComponent], toleranceV));
#endif
                } // for
            } // for
        } // for
    } // for

    PYLITH_METHOD_END;
} // _checkVertexFields


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_TestDataWriterHDF5::MeshData::_checkCellFields(const MeshData& reference,
                                                                const MeshData& candidate) {
    PYLITH_METHOD_BEGIN;

#if 0 // TEMPORARY
    REQUIRE(candidate.vertexFields.size() == reference.vertexFields.size());
#endif
    for (const auto& kv : reference.vertexFields) {
        const std::string& name = kv.first;
        const Field& referenceField = kv.second;

        const auto& iterC = candidate.vertexFields.find(name);
        REQUIRE(iterC != candidate.vertexFields.end());
        const Field& candidateField = iterC->second;

        INFO("Checking cell field '" << name << "'.");

        REQUIRE(reference.numVertices == candidate.numVertices);
        REQUIRE(reference.spaceDim == candidate.spaceDim);
        REQUIRE(reference.numCells == candidate.numCells);
        REQUIRE(reference.numCorners == candidate.numCorners);
        REQUIRE(referenceField.numPoints == candidateField.numPoints);
        REQUIRE(referenceField.numComponents == candidateField.numComponents);

        const size_t numComponents = referenceField.numComponents;
        const size_t spaceDim = reference.spaceDim;
        const size_t numCorners = reference.numCorners;
        const std::vector<size_t>& referenceOrder = reference._getCellsOrder();
        const std::vector<size_t>& candidateOrder = candidate._getCellsOrder();

        const PylithScalar tolerance = 1.0e-6;
        for (size_t iCell = 0; iCell < reference.numCells; ++iCell) {
            const size_t cellReference = referenceOrder[iCell];
            const size_t cellCandidate = candidateOrder[iCell];
            INFO("Checking candidate cell " << cellCandidate << " against reference cell " << cellReference << ".");

            for (size_t iComponent = 0; iComponent < numComponents; ++iComponent) {
                const PylithScalar toleranceV = std::max(tolerance, tolerance*referenceField.data[cellReference*numComponents+iComponent]);
#if 0
                CHECK_THAT(candidateField.data[cellCandidate*numComponents+iComponent],
                           Catch::Matchers::WithinAbs(referenceField.data[cellReference*numComponents+iComponent], toleranceV));
#endif
            } // for
        } // for
    } // for

    PYLITH_METHOD_END;
} // _checkCellFields


// ------------------------------------------------------------------------------------------------
herr_t
pylith_meshio_TestDataWriterHDF5_checkObject(hid_t id,
                                             const char* name,
                                             const H5O_info_t* info,
                                             void* data) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(info);
    REQUIRE(data);

    INFO("Checking dataset '" << name << "'.");

    hid_t* file = (hid_t*) data;REQUIRE(H5Iis_valid(*file));
    herr_t err = 0;

    switch (info->type) {
    case H5O_TYPE_GROUP: {
        hid_t group = H5Gopen2(*file, name, H5P_DEFAULT);REQUIRE(group >= 0);
        err = H5Gclose(group);REQUIRE(err >= 0);
        break;
    } // group
    case H5O_TYPE_DATASET: {
        // Get expected dataset.
        hid_t datasetE = H5Dopen2(id, name, H5P_DEFAULT);REQUIRE(datasetE >= 0);
        hid_t dataspaceE = H5Dget_space(datasetE);REQUIRE(dataspaceE >= 0);
        const int ndimsE = H5Sget_simple_extent_ndims(dataspaceE);REQUIRE(ndimsE > 0);
        hsize_t* dimsE = (ndimsE > 0) ? new hsize_t[ndimsE] : 0;
        const int ndimsECheck = H5Sget_simple_extent_dims(dataspaceE, dimsE, 0);
        REQUIRE(ndimsE == ndimsECheck);
        int sizeE = (ndimsE > 0 && dimsE[0] > 0) ? 1 : 0;
        for (int i = 0; i < ndimsE; ++i) {
            sizeE *= dimsE[i];
        } // for

        // Get dataset
        hid_t dataset = H5Dopen2(*file, name, H5P_DEFAULT);REQUIRE(dataset >= 0);
        hid_t dataspace = H5Dget_space(dataset);REQUIRE(dataspace >= 0);
        const int ndims = H5Sget_simple_extent_ndims(dataspace);REQUIRE(ndims > 0);
        hsize_t* dims = (ndims > 0) ? new hsize_t[ndims] : 0;
        const int ndimsCheck = H5Sget_simple_extent_dims(dataspace, dims, 0);
        REQUIRE(ndims == ndimsCheck);
        int size = (ndims > 0 && dims[0] > 0) ? 1 : 0;
        for (int i = 0; i < ndims; ++i) {
            size *= dims[i];
        } // for

        // Check dimensions.
        REQUIRE(ndimsE == ndims);
        for (int i = 0; i < ndimsE; ++i) {
            CHECK(dimsE[i] == dims[i]);
        } // for

        // Check the expected datatype
        hid_t datatypeE = H5Dget_type(datasetE);REQUIRE(datatypeE >= 0);
        hid_t dataclassE = H5Tget_class(datatypeE);REQUIRE(dataclassE >= 0);

        hid_t datatype = H5Dget_type(dataset);REQUIRE(datatype >= 0);
        hid_t dataclass = H5Tget_class(datatype);REQUIRE(dataclass >= 0);

        switch (dataclassE) {
        case H5T_FLOAT: {
            double* dataE = (sizeE > 0) ? new double[sizeE] : 0;REQUIRE(sizeE > 0);
            err = H5Dread(datasetE, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*) dataE);REQUIRE(err >= 0);

            double* data = (size > 0) ? new double[size] : 0;REQUIRE(size > 0);
            err = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*) data);REQUIRE(err >= 0);

            REQUIRE(sizeE == size);

            // Compare data values.
            const double tolerance = 1.0e-6;
            for (int i = 0; i < size; ++i) {
                const double toleranceV = std::max(tolerance, tolerance*dataE[i]);
                CHECK_THAT(double(data[i]), Catch::Matchers::WithinAbs(dataE[i], toleranceV));
            } // for

            delete[] dataE;dataE = 0;
            delete[] data;data = 0;

            break;
        } // H5T_DOUBLE

        case H5T_INTEGER: {
            int* dataE = (sizeE > 0) ? new int[sizeE] : 0;REQUIRE(sizeE > 0);
            err = H5Dread(datasetE, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*) dataE);REQUIRE(err >= 0);

            int* data = (size > 0) ? new int[size] : 0;REQUIRE(size > 0);
            err = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*) data);REQUIRE(err >= 0);

            REQUIRE(sizeE == size);

            // Compare data values.
            for (int i = 0; i < size; ++i) {
                CHECK(dataE[i] == int(data[i]));
            } // for

            delete[] dataE;dataE = 0;
            delete[] data;data = 0;

            break;
        } // H5T_INTEGER

        case H5T_STRING: {
            const int slenE = H5Tget_size(datatypeE);REQUIRE(slenE > 0);
            sizeE *= slenE;

            const int slen = H5Tget_size(datatype);REQUIRE(slen > 0);
            size *= slen;

            REQUIRE(slenE == slen);
            REQUIRE(sizeE == size);

            char* dataE = (sizeE > 0) ? new char[sizeE] : 0;REQUIRE(sizeE > 0);
            err = H5Dread(datasetE, datatypeE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataE);REQUIRE(err >= 0);

            char* data = (size > 0) ? new char[size] : 0;REQUIRE(size > 0);
            err = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);REQUIRE(err >= 0);

            for (int i = 0; i < size; ++i) {
                CHECK(dataE[i] == data[i]);
            } // for

            delete[] dataE;dataE = 0;
            delete[] data;data = 0;

            break;
        } // H5T_C_S1

        default:
            REQUIRE(false);
        } // switch

        err = H5Sclose(dataspaceE);REQUIRE(err >= 0);
        err = H5Dclose(datasetE);REQUIRE(err >= 0);

        err = H5Sclose(dataspace);REQUIRE(err >= 0);
        err = H5Dclose(dataset);REQUIRE(err >= 0);

        delete[] dimsE;dimsE = 0;
        delete[] dims;dims = 0;

        break;
    } // dataset
    default:
        REQUIRE(false);
    } // switch

    PYLITH_METHOD_RETURN(0);
} // checkObject


// ------------------------------------------------------------------------------------------------
// Check HDF5 file against archived file.
void
pylith::meshio::TestDataWriterHDF5::checkFile(const char* filename) {
    PYLITH_METHOD_BEGIN;

    const std::string filenameE = "data/" + std::string(filename);

#if 1
    std::unique_ptr<_TestDataWriterHDF5::MeshData> reference(_TestDataWriterHDF5::MeshData::load(filenameE.c_str()));
    std::unique_ptr<_TestDataWriterHDF5::MeshData> candidate(_TestDataWriterHDF5::MeshData::load(filename));

    _TestDataWriterHDF5::MeshData::verify(*reference, *candidate);
#else
    herr_t err = 0;

    hid_t fileE = H5Fopen(filenameE.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);REQUIRE(fileE >= 0);

    hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);REQUIRE(file >= 0);
#if defined(PYLITH_HDF5_USE_API_112)
    // Traverse recursively file with expected values.
    err = H5Ovisit(fileE, H5_INDEX_NAME, H5_ITER_NATIVE, pylith_meshio_TestDataWriterHDF5_checkObject, (void*) &file, H5O_INFO_ALL);REQUIRE(err >= 0);
#else
    err = H5Ovisit(fileE, H5_INDEX_NAME, H5_ITER_NATIVE, pylith_meshio_TestDataWriterHDF5_checkObject, (void*) &file);REQUIRE(err >= 0);
#endif
    err = H5Fclose(fileE);
    REQUIRE(err >= 0);

    err = H5Fclose(file);
    REQUIRE(err >= 0);
#endif


    PYLITH_METHOD_END;
} // checkFile


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::TestDataWriterHDF5_Data::TestDataWriterHDF5_Data(void) :
    opencloseFilename(NULL),
    vertexFilename(NULL),
    cellFilename(NULL) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::TestDataWriterHDF5_Data::~TestDataWriterHDF5_Data(void) {}


// End of file
