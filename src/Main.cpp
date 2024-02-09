/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2023-2024, Christoph Neuhauser
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include <unordered_map>
#include <filesystem>

#include <glm/glm.hpp>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Unstable/Topology/TetTopology.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>

#include "LineReader.hpp"
#include "CSPSolver.hpp"
#include "NaxosSolver.hpp"
#include "FlipSolver.hpp"

/*
 * b: 0, c: 1, d: 2, b': 3, c': 4, d': 5
 * abc: R(0) -> (b, c'), F(1) -> (c, b')
 * acd: R(0) -> (c, d'), F(1) -> (d, c')
 * abd: R(0) -> (d, b'), F(1) -> (b, d')
 */
const uint32_t PRISM_TO_TET_TABLE[8][12] = {
        // RRR, invalid
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        // RRF, splits: (b, c'), (c, d'), (b, d'); meets in: b, d'
        { 0, 1, 5, 2, 3, 4, 5, 0, 0, 1, 4, 5 },
        // RFR, splits: (b, c'), (d, c'), (d, b'); meets in: d, c'
        { 0, 1, 4, 2, 3, 4, 5, 2, 0, 2, 4, 3 },
        // RFF, splits: (b, c'), (d, c'), (b, d'); meets in: b, c'
        { 0, 1, 4, 2, 3, 4, 5, 0, 0, 2, 4, 5 },
        // FRR, splits: (c, b'), (c, d'), (d, b'); meets in: c, b'
        { 0, 1, 3, 2, 3, 4, 5, 1, 1, 2, 5, 3 },
        // FRF, splits: (c, b'), (c, d'), (b, d'); meets in: c, d'
        { 0, 1, 5, 2, 3, 4, 5, 1, 0, 1, 3, 5 },
        // FFR, splits: (c, b'), (d, c'), (d, b'); meets in: d, b'
        { 0, 1, 3, 2, 3, 4, 5, 2, 1, 2, 4, 3 },
        // FFF, invalid
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
// Original before swaps for winding:
/*const uint32_t PRISM_TO_TET_TABLE[8][12] = {
        // RRR, invalid
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        // RRF, splits: (b, c'), (c, d'), (b, d'); meets in: b, d'
        { 0, 1, 2, 5, 3, 4, 5, 0, 0, 1, 4, 5 },
        // RFR, splits: (b, c'), (d, c'), (d, b'); meets in: d, c'
        { 0, 1, 2, 4, 3, 4, 5, 2, 0, 2, 3, 4 },
        // RFF, splits: (b, c'), (d, c'), (b, d'); meets in: b, c'
        { 0, 1, 2, 4, 3, 4, 5, 0, 0, 2, 4, 5 },
        // FRR, splits: (c, b'), (c, d'), (d, b'); meets in: c, b'
        { 0, 1, 2, 3, 3, 4, 5, 1, 1, 2, 3, 5 },
        // FRF, splits: (c, b'), (c, d'), (b, d'); meets in: c, d'
        { 0, 1, 2, 5, 3, 4, 5, 1, 0, 1, 3, 5 },
        // FFR, splits: (c, b'), (d, c'), (d, b'); meets in: d, b'
        { 0, 1, 2, 3, 3, 4, 5, 2, 1, 2, 3, 4 },
        // FFF, invalid
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};*/

class TetMesh {
private:
    //OpenVolumeMesh::GeometricTetrahedralMeshV3f ovmMesh;
    OpenVolumeMesh::TetrahedralGeometryKernel<OpenVolumeMesh::Geometry::Vec3f, OpenVolumeMesh::TetrahedralMeshTopologyKernel> ovmMesh;

public:
    size_t getNumFaces() {
        return ovmMesh.n_faces();
    }

    // Build a tetrahedral mesh from a list of cell indices (4 vertex indices define a tet) and vertex positions.
    void build(const std::vector<uint32_t>& cellIndices, const std::vector<glm::vec3>& vertexPositions) {
        // https://www.graphics.rwth-aachen.de/media/openvolumemesh_static/Documentation/OpenVolumeMesh-Doc-Latest/ovm_tutorial_01.html
        ovmMesh = {};
        std::vector<OpenVolumeMesh::VertexHandle> ovmVertices;
        ovmVertices.reserve(vertexPositions.size());
        for (const glm::vec3& v : vertexPositions) {
            ovmVertices.emplace_back(ovmMesh.add_vertex(OpenVolumeMesh::Vec3f(v.x, v.y, v.z)));
        }

        const uint32_t numTets = cellIndices.size() / 4;
        std::vector<OpenVolumeMesh::VertexHandle> ovmCellVertices(4);
        for (uint32_t cellId = 0; cellId < numTets; ++cellId) {
            for (uint32_t faceIdx = 0; faceIdx < 4; faceIdx++){
                ovmCellVertices.at(faceIdx) = ovmVertices.at(cellIndices.at(cellId * 4 + faceIdx));
            }
            ovmMesh.add_cell(ovmCellVertices);
        }
    }

    bool loadTxt(const std::string& filePath) {
        std::vector<uint32_t> cellIndices;
        std::vector<glm::vec3> vertexPositions;

        uint8_t* buffer = nullptr;
        size_t length = 0;
        bool loaded = loadFileFromSource(filePath, buffer, length, true);
        if (!loaded) {
            return false;
        }

        LineReader lineReader(reinterpret_cast<const char*>(buffer), length);
        std::vector<std::string> linesInfo;
        std::vector<uint32_t> cellIndicesVector;
        std::vector<float> vertexPositionsVector;
        std::vector<float> vertexColorsVector;
        while (lineReader.isLineLeft()) {
            lineReader.readVectorLine<std::string>(linesInfo);
            if (linesInfo.size() != 2) {
                throw std::runtime_error("Error in TxtTetLoader::loadFromFile: Invalid header line.");
                delete[] buffer;
                return false;
            }
            const std::string& key = linesInfo.at(0);
            const auto numEntries = fromString<uint32_t>(linesInfo.at(1));
            if (key == "cellIndices") {
                cellIndices.reserve(4 * numEntries);
                for (uint32_t i = 0; i < numEntries; i++) {
                    lineReader.readVectorLine<uint32_t>(cellIndicesVector);
                    if (cellIndicesVector.size() != 4) {
                        throw std::runtime_error("Error in TxtTetLoader::loadFromFile: Invalid number of indices.");
                        delete[] buffer;
                        return false;
                    }
                    for (auto entry : cellIndicesVector) {
                        cellIndices.push_back(entry);
                    }
                }
            } else if (key == "vertexPositions") {
                vertexPositions.reserve(numEntries);
                for (uint32_t i = 0; i < numEntries; i++) {
                    lineReader.readVectorLine<float>(vertexPositionsVector);
                    if (vertexPositionsVector.size() != 3) {
                        throw std::runtime_error(
                                "Error in TxtTetLoader::loadFromFile: Invalid number of vertex position entries.");
                        delete[] buffer;
                        return false;
                    }
                    vertexPositions.emplace_back(
                            vertexPositionsVector.at(0), vertexPositionsVector.at(1), vertexPositionsVector.at(2));
                }
            } else if (key == "vertexColors") {
                //vertexColors.reserve(numEntries);
                for (uint32_t i = 0; i < numEntries; i++) {
                    lineReader.readVectorLine<float>(vertexColorsVector);
                    if (vertexColorsVector.size() != 4) {
                        throw std::runtime_error(
                                "Error in TxtTetLoader::loadFromFile: Invalid number of vertex color entries.");
                        delete[] buffer;
                        return false;
                    }
                    // Skip.
                    //vertexColors.emplace_back(
                    //        vertexColorsVector.at(0), vertexColorsVector.at(1),
                    //        vertexColorsVector.at(2), vertexColorsVector.at(3));
                }
            }
        }

        delete[] buffer;

        build(cellIndices, vertexPositions);
        return true;
    }

#define USE_SPLIT_EDGE

    /**
     * Subdivides all tetrahedra incident with the specified vertex.
     * The edges incident with the vertex are split by factor t \in (0, 1).
     * The old cells are removed, and new cells are added.
     */
    void subdivideAtVertex(uint32_t vertexIndex, float t) {
        auto& vertexPositionsOvm = ovmMesh.vertex_positions();

        // Add new vertices along the incident edges & collect old edges to delete.
        std::vector<OpenVolumeMesh::EdgeHandle> edgesToDelete;
        std::unordered_map<OpenVolumeMesh::EdgeHandle, OpenVolumeMesh::VertexHandle> edgeToNewVertexMap;
        OpenVolumeMesh::VertexHandle vh((int)vertexIndex);
        for (auto ve_it = ovmMesh.ve_iter(vh); ve_it.valid(); ve_it++) {
            const auto& eh = ve_it.cur_handle();
            edgesToDelete.push_back(eh);
#ifndef USE_SPLIT_EDGE
            auto vhs = ovmMesh.edge_vertices(eh);
            auto vh0 = vhs.at(0);
            auto vh1 = vhs.at(1);
            const auto& vp0Ovm = vertexPositionsOvm.at(vh0);
            const auto& vp1Ovm = vertexPositionsOvm.at(vh1);
            glm::vec3 vp0(vp0Ovm[0], vp0Ovm[1], vp0Ovm[2]);
            glm::vec3 vp1(vp1Ovm[0], vp1Ovm[1], vp1Ovm[2]);
            glm::vec3 vpe = glm::mix(vp0, vp1, 0.5f);
            auto vhe = ovmMesh.add_vertex(OpenVolumeMesh::Vec3f(vpe.x, vpe.y, vpe.z));
            edgeToNewVertexMap.insert(std::make_pair(eh, vhe));
#endif
        }

#ifdef USE_SPLIT_EDGE
        for (const auto& eh : edgesToDelete) {
            ovmMesh.split_edge(eh);
        }
        ovmMesh.collect_garbage();
#else
        /*
         * Collect the new indices of the new cells. Only add them after deleting the old cells.
         *
         * Should adding the new cells be OK before deleting the old cells?
         * The mesh would not be manifold in the intermediate state.
         */
        std::vector<OpenVolumeMesh::VertexHandle> newCells;

        // Each cell has 3 subdivided edges. 'vhes' stores the 3 subdivision points.
        std::array<OpenVolumeMesh::VertexHandle, 3> vhes;
        // 'vhbs' stores the end points of the 3 subdivided edges not equal to 'vh'.
        std::array<OpenVolumeMesh::VertexHandle, 3> vhbs;
        // Iterate over all cells incident with the vertex.
        for (auto vc_it = ovmMesh.vc_iter(vh); vc_it.valid(); vc_it++) {
            //auto incidentVertices = ovmMesh.get_cell_vertices(*vc_it, vh);

            auto tt = OpenVolumeMesh::TetTopology(ovmMesh, *vc_it, vh); // vertex a == vh
            vhes[0] = edgeToNewVertexMap.find(tt.ab().edge_handle())->second;
            vhbs[0] = tt.b();
            vhes[1] = edgeToNewVertexMap.find(tt.ac().edge_handle())->second;
            vhbs[1] = tt.c();
            vhes[2] = edgeToNewVertexMap.find(tt.ad().edge_handle())->second;
            vhbs[2] = tt.d();

            // Iterate over all edges incident with the cell.
            /*int i = 0;
            for (auto ce_it = ovmMesh.ce_iter(vc_it.cur_handle()); ce_it.valid(); ce_it++) {
                auto vhs = ovmMesh.edge_vertices(ce_it.cur_handle());
                // Only process this edge if it was subdivided (i.e., is incident with the vertex 'vh').
                auto ev_it = edgeToNewVertexMap.find(ce_it.cur_handle());
                if (ev_it == edgeToNewVertexMap.end()) {
                    continue;
                }
                OpenVolumeMesh::VertexHandle vh0 = vhs[0];
                OpenVolumeMesh::VertexHandle vh1 = vhs[1];
                // Retrieve the new subdivision vertex & the end point of the subdivided edge not equal to 'vh'.
                vhes[i] = ev_it->second;
                vhbs[i] = vh == vh0 ? vh1 : vh0;
                i++;
            }*/

            // Collect the new cell vertex indices subdividing the currently iterated cell.
            newCells.push_back(vh);
            newCells.push_back(vhes[0]);
            newCells.push_back(vhes[1]);
            newCells.push_back(vhes[2]);

            newCells.push_back(vhes[0]);
            newCells.push_back(vhbs[0]);
            newCells.push_back(vhbs[1]);
            newCells.push_back(vhbs[2]);

            newCells.push_back(vhes[0]);
            newCells.push_back(vhes[1]);
            newCells.push_back(vhbs[2]);
            newCells.push_back(vhbs[1]);

            newCells.push_back(vhes[0]);
            newCells.push_back(vhes[1]);
            newCells.push_back(vhes[2]);
            newCells.push_back(vhbs[2]);
        }

        // If we sort edge handles from largest to smallest, there should hopefully not be a problem with invalid handles.
        std::sort(edgesToDelete.rbegin(), edgesToDelete.rend());
        for (const auto& eh : edgesToDelete) {
            ovmMesh.delete_edge(eh);
        }
        ovmMesh.collect_garbage();

        // Add the new cells after the old ones have been deleted.
        for (size_t i = 0; i < newCells.size(); i += 4) {
            auto& p0Ovm = vertexPositionsOvm.at(newCells.at(i));
            auto& p1Ovm = vertexPositionsOvm.at(newCells.at(i + 1));
            auto& p2Ovm = vertexPositionsOvm.at(newCells.at(i + 2));
            auto& p3Ovm = vertexPositionsOvm.at(newCells.at(i + 3));
            glm::vec3 p0(p0Ovm[0], p0Ovm[1], p0Ovm[2]);
            glm::vec3 p1(p1Ovm[0], p1Ovm[1], p1Ovm[2]);
            glm::vec3 p2(p2Ovm[0], p2Ovm[1], p2Ovm[2]);
            glm::vec3 p3(p3Ovm[0], p3Ovm[1], p3Ovm[2]);
            float volumeSign = -glm::sign(glm::dot(glm::cross(p1 - p0, p2 - p0), p3 - p0));
            assert(volumeSign > 0.0f && "Invalid winding");
            ovmMesh.add_cell(newCells.at(i), newCells.at(i + 1), newCells.at(i + 2), newCells.at(i + 3), true);
        }
#endif

        // Sanity check: Every tetrahedral cell may only have 4 vertices.
        for (OpenVolumeMesh::CellIter c_it = ovmMesh.cells_begin(); c_it != ovmMesh.cells_end(); c_it++) {
            auto ch = *c_it;
            int vidx = 0;
            for (auto cv_it = ovmMesh.cv_iter(ch); cv_it.valid(); cv_it++) {
                if (vidx >= 4) {
                    int i = 0;
                    for (auto cv_it = ovmMesh.cv_iter(ch); cv_it.valid(); cv_it++) {
                        std::cout << "v" << i << ": " << cv_it->uidx() << std::endl;
                        i++;
                    }
                }
                assert(vidx < 4);
                vidx++;
            }
        }

        // Sanity check: Winding.
        for (OpenVolumeMesh::CellIter c_it = ovmMesh.cells_begin(); c_it != ovmMesh.cells_end(); c_it++) {
            auto ch = *c_it;
            auto cellVertices = ovmMesh.get_cell_vertices(ch);
            auto& p0Ovm = vertexPositionsOvm.at(cellVertices.at(0));
            auto& p1Ovm = vertexPositionsOvm.at(cellVertices.at(1));
            auto& p2Ovm = vertexPositionsOvm.at(cellVertices.at(2));
            auto& p3Ovm = vertexPositionsOvm.at(cellVertices.at(3));
            glm::vec3 p0(p0Ovm[0], p0Ovm[1], p0Ovm[2]);
            glm::vec3 p1(p1Ovm[0], p1Ovm[1], p1Ovm[2]);
            glm::vec3 p2(p2Ovm[0], p2Ovm[1], p2Ovm[2]);
            glm::vec3 p3(p3Ovm[0], p3Ovm[1], p3Ovm[2]);
            float volumeSign = -glm::sign(glm::dot(glm::cross(p1 - p0, p2 - p0), p3 - p0));
            assert(volumeSign > 0.0f && "Invalid winding");
        }
    }


    void subdivideAtVertexPrism(uint32_t vertexIndex, float t) {
        auto& vertexPositionsOvm = ovmMesh.vertex_positions();
        OpenVolumeMesh::VertexHandle vh((int)vertexIndex);

        // Add new vertices along the incident edges & collect old edges to delete.
        std::vector<OpenVolumeMesh::EdgeHandle> edgesToDelete;
        std::unordered_map<OpenVolumeMesh::EdgeHandle, OpenVolumeMesh::VertexHandle> edgeToNewVertexMap;
        for (auto ve_it = ovmMesh.ve_iter(vh); ve_it.valid(); ve_it++) {
            const auto& eh = ve_it.cur_handle();
            edgesToDelete.push_back(eh);
            auto vhs = ovmMesh.edge_vertices(eh);
            auto vh0 = vhs.at(0);
            auto vh1 = vhs.at(1);
            const auto& vp0Ovm = vertexPositionsOvm.at(vh0);
            const auto& vp1Ovm = vertexPositionsOvm.at(vh1);
            glm::vec3 vp0(vp0Ovm[0], vp0Ovm[1], vp0Ovm[2]);
            glm::vec3 vp1(vp1Ovm[0], vp1Ovm[1], vp1Ovm[2]);
            glm::vec3 vpe = glm::mix(vp0, vp1, 0.5f);
            auto vhe = ovmMesh.add_vertex(OpenVolumeMesh::Vec3f(vpe.x, vpe.y, vpe.z));
            edgeToNewVertexMap.insert(std::make_pair(eh, vhe));
        }

        // Iterate over all cells incident with the vertex.
        int numIncidentTets = 0;
        for (auto vc_it = ovmMesh.vc_iter(vh); vc_it.valid(); vc_it++) {
            numIncidentTets++;
        }
        auto** tetTopologies = new OpenVolumeMesh::TetTopology*[numIncidentTets];
        int tetIdx = 0;
        std::unordered_map<OpenVolumeMesh::HFH, int> halffaceToIndexMap;
        std::unordered_map<OpenVolumeMesh::CH, int> cellToIndexMap;
        for (auto vc_it = ovmMesh.vc_iter(vh); vc_it.valid(); vc_it++) {
            auto* tt = new OpenVolumeMesh::TetTopology(ovmMesh, *vc_it, vh); // vertex a == vh
            tetTopologies[tetIdx] = tt;
            halffaceToIndexMap.insert(std::make_pair(tt->abc(), tetIdx * 3));
            halffaceToIndexMap.insert(std::make_pair(tt->acd(), tetIdx * 3 + 1));
            halffaceToIndexMap.insert(std::make_pair(tt->adb(), tetIdx * 3 + 2));
            cellToIndexMap.insert(std::make_pair(*vc_it, tetIdx));
            tetIdx++;
        }

        // Build neighborhood graph.
        std::vector<Prism> prisms(numIncidentTets);
        for (tetIdx = 0; tetIdx < numIncidentTets; tetIdx++) {
            auto& prism = prisms.at(tetIdx);
            auto& tt = *tetTopologies[tetIdx]; // vertex a == vh
            std::array<OpenVolumeMesh::HFH, 3> hfhs = {
                    tt.abc().opposite_handle(), tt.acd().opposite_handle(), tt.adb().opposite_handle()
            };
            for (int faceIdx = 0; faceIdx < 3; faceIdx++) {
                if (hfhs.at(faceIdx).is_valid()) {
                    auto itFace = halffaceToIndexMap.find(hfhs.at(faceIdx));
                    auto ch = ovmMesh.incident_cell(hfhs.at(faceIdx));
                    auto itCell = cellToIndexMap.find(ch);
                    prism.neighborFaceIndices.at(faceIdx) = itFace != halffaceToIndexMap.end() ? itFace->second : -1;
                    prism.neighbors.at(faceIdx) = itCell != cellToIndexMap.end() ? itCell->second : -1;
                }
            }
        }

        // Solve constraint satisfaction problem.
        //auto* cspSolver = new NaxosSolver();
        auto* cspSolver = new FlipSolver();
        if (!cspSolver->solve(prisms)) {
            throw std::runtime_error("Error: CSP solver failed!");
        }
        if (!checkIsCspFulfilled(prisms)) {
            throw std::runtime_error("Error: CSP condition not fulfilled!");
        }
        delete cspSolver;

        // Collect the new indices of the new cells. Only add them after deleting the old cells.
        std::vector<OpenVolumeMesh::VertexHandle> newCells;
        // Prism vertex handles.
        std::array<OpenVolumeMesh::VertexHandle, 6> pvhs;
        for (tetIdx = 0; tetIdx < numIncidentTets; tetIdx++) {
            auto& prism = prisms.at(tetIdx);
            auto& tt = *tetTopologies[tetIdx]; // vertex a == vh
            pvhs[0] = tt.b();
            pvhs[1] = tt.c();
            pvhs[2] = tt.d();
            pvhs[3] = edgeToNewVertexMap.find(tt.ab().edge_handle())->second;
            pvhs[4] = edgeToNewVertexMap.find(tt.ac().edge_handle())->second;
            pvhs[5] = edgeToNewVertexMap.find(tt.ad().edge_handle())->second;

            // Add top.
            newCells.push_back(vh);
            newCells.push_back(pvhs[3]);
            newCells.push_back(pvhs[4]);
            newCells.push_back(pvhs[5]);

            uint32_t splits = prism.cuts.bitfield;
            for (int i = 0; i < 12; i++) {
                newCells.push_back(pvhs[PRISM_TO_TET_TABLE[splits][i]]);
            }
        }

        for (tetIdx = 0; tetIdx < numIncidentTets; tetIdx++) {
            delete tetTopologies[tetIdx];
        }
        delete[] tetTopologies;

        for (const auto& eh : edgesToDelete) {
            ovmMesh.delete_edge(eh);
        }
        ovmMesh.collect_garbage();

        // Add the new cells after the old ones have been deleted.
        for (size_t i = 0; i < newCells.size(); i += 4) {
            auto& p0Ovm = vertexPositionsOvm.at(newCells.at(i));
            auto& p1Ovm = vertexPositionsOvm.at(newCells.at(i + 1));
            auto& p2Ovm = vertexPositionsOvm.at(newCells.at(i + 2));
            auto& p3Ovm = vertexPositionsOvm.at(newCells.at(i + 3));
            glm::vec3 p0(p0Ovm[0], p0Ovm[1], p0Ovm[2]);
            glm::vec3 p1(p1Ovm[0], p1Ovm[1], p1Ovm[2]);
            glm::vec3 p2(p2Ovm[0], p2Ovm[1], p2Ovm[2]);
            glm::vec3 p3(p3Ovm[0], p3Ovm[1], p3Ovm[2]);
            float volumeSign = -glm::sign(glm::dot(glm::cross(p1 - p0, p2 - p0), p3 - p0));
            assert(volumeSign > 0.0f && "Invalid winding");
            ovmMesh.add_cell(newCells.at(i), newCells.at(i + 1), newCells.at(i + 2), newCells.at(i + 3), true);
        }
    }

    void checkPrismTableWinding(uint32_t tetIdx, uint32_t vertexIndex, float t = 0.5f) {
        auto& vertexPositionsOvm = ovmMesh.vertex_positions();
        OpenVolumeMesh::VertexHandle vh((int)vertexIndex);
        OpenVolumeMesh::CellHandle ch((int)tetIdx);

        OpenVolumeMesh::TetTopology tt(ovmMesh, ch, vh); // vertex a == vh
        const auto& vpaOvm = vertexPositionsOvm.at(tt.a());
        const auto& vpbOvm = vertexPositionsOvm.at(tt.b());
        const auto& vpcOvm = vertexPositionsOvm.at(tt.c());
        const auto& vpdOvm = vertexPositionsOvm.at(tt.d());
        glm::vec3 vpa(vpaOvm[0], vpaOvm[1], vpaOvm[2]);
        glm::vec3 vpb(vpbOvm[0], vpbOvm[1], vpbOvm[2]);
        glm::vec3 vpc(vpcOvm[0], vpcOvm[1], vpcOvm[2]);
        glm::vec3 vpd(vpdOvm[0], vpdOvm[1], vpdOvm[2]);
        glm::vec3 vpbp = glm::mix(vpa, vpb, 0.5f);
        glm::vec3 vpcp = glm::mix(vpa, vpc, 0.5f);
        glm::vec3 vpdp = glm::mix(vpa, vpd, 0.5f);

        std::array<glm::vec3, 6> prismVertexPositions{};
        prismVertexPositions[0] = vpb;
        prismVertexPositions[1] = vpc;
        prismVertexPositions[2] = vpd;
        prismVertexPositions[3] = vpbp;
        prismVertexPositions[4] = vpcp;
        prismVertexPositions[5] = vpdp;

        for (int caseIdx = 1; caseIdx < 7; caseIdx++) {
            for (int tetIdxLocal = 0; tetIdxLocal < 3; tetIdxLocal++) {
                glm::vec3 p0 = prismVertexPositions[PRISM_TO_TET_TABLE[caseIdx][tetIdxLocal * 4 + 0]];
                glm::vec3 p1 = prismVertexPositions[PRISM_TO_TET_TABLE[caseIdx][tetIdxLocal * 4 + 1]];
                glm::vec3 p2 = prismVertexPositions[PRISM_TO_TET_TABLE[caseIdx][tetIdxLocal * 4 + 2]];
                glm::vec3 p3 = prismVertexPositions[PRISM_TO_TET_TABLE[caseIdx][tetIdxLocal * 4 + 3]];
                float signVal = glm::sign(glm::dot(glm::cross(p1 - p0, p2 - p0), p3 - p0));
                if (signVal >= 0.0f) {
                    std::cout << "Invalid sign for case " << caseIdx << ", tet " << tetIdxLocal << std::endl;
                }
            }
        }
    }
};

int main() {
    TetMesh tetMesh;
    std::string dataDir = "Data/";
    if (!std::filesystem::is_directory(dataDir)) {
        dataDir = "../Data/";
    }

    int testCaseIdx = 1;
    if (testCaseIdx == 0) {
        // Two tetrahedra sharing one face.
        std::vector<uint32_t> cellIndices = {
                0, 1, 2, 3,
                1, 4, 2, 3,
        };
        std::vector<glm::vec3> vertexPositions = {
                { -0.1, -0.1,  0.1 },
                {  0.1, -0.1,  0.1 },
                {    0,  0.3,    0 },
                {  0.1, -0.1, -0.1 },
                {  0.2,  0.1,  0.1 },
        };
        tetMesh.build(cellIndices, vertexPositions);
        for (int vertexIdx = 0; vertexIdx < 3; vertexIdx++) {
            tetMesh.checkPrismTableWinding(0, vertexIdx);
        }
        for (int vertexIdx = 1; vertexIdx < 4; vertexIdx++) {
            tetMesh.checkPrismTableWinding(1, vertexIdx);
        }
        tetMesh.subdivideAtVertexPrism(0, 0.5f);
    } else if (testCaseIdx == 1) {
        tetMesh.loadTxt(dataDir + "ico01.txt");
        tetMesh.subdivideAtVertexPrism(12, 0.5f);
        //tetMesh.subdivideAtVertexPrism(1, 0.5f);
    }

//#ifdef USE_SPLIT_EDGE
//    std::cout << "#Faces before: " << tetMesh.getNumFaces() << std::endl;
//    tetMesh.subdivideAtVertex(4, 0.5f);
//    std::cout << "#Faces after: " << tetMesh.getNumFaces() << std::endl;
//#else
//    for (int j = 0; j < 4; j++) {
//        tetMesh.subdivideAtVertex(j, 0.5f);
//    }
//#endif

    return 0;
}
