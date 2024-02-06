/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2024, Christoph Neuhauser
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

#include <naxos.h>
#include "NaxosSolver.hpp"

bool NaxosSolver::solve(std::vector<Prism>& prisms) {
    naxos::NsProblemManager pm;
    naxos::NsIntVarArray vars;
    //std::vector<naxos::NsIntVar*> varsVector(prisms.size() * 3);

    // Add all variables.
    for (int prismIdx = 0; prismIdx < int(prisms.size()); prismIdx++) {
        for (int faceIdx = 0; faceIdx < 3; faceIdx++) {
            vars.push_back(naxos::NsIntVar(pm, 0, 1));
            //varsVector.at(3 * prismIdx + faceIdx) = &vars.back();
        }
    }

    // Add the constraint that no tet may only have rising or falling cuts.
    for (int prismIdx = 0; prismIdx < int(prisms.size()); prismIdx++) {
        naxos::NsIntVar sum = naxos::NsSum(vars, prismIdx * 3, 3);
        pm.add(sum != 0);
        pm.add(sum != 3);
    }

    // Add the constraint that neighboring tets must have opposite cuts.
    for (int prismIdx = 0; prismIdx < int(prisms.size()); prismIdx++) {
        auto& prism = prisms.at(prismIdx);
        for (int faceIdx = 0; faceIdx < 3; faceIdx++) {
            int neighborFaceIdx = prism.neighborFaceIndices.at(faceIdx);
            if (neighborFaceIdx != -1) {
                //pm.add(*varsVector.at(prismIdx * 3 + faceIdx) != *varsVector.at(neighborFaceIdx));
            }
        }
    }

    // Solve the constraint satisfaction problem.
    pm.addGoal(new naxos::NsgLabeling(vars));
    if (pm.nextSolution()) {
        for (int prismIdx = 0; prismIdx < int(prisms.size()); prismIdx++) {
            auto& prism = prisms.at(prismIdx);
            for (int faceIdx = 0; faceIdx < 3; faceIdx++) {
                //prism.cuts.setCut(uint32_t(faceIdx), uint32_t(varsVector.at(3 * prismIdx + faceIdx)->value()));
            }
        }
        return true;
    }
    return false;
}
