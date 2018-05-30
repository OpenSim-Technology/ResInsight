/////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (C) 2018-     Statoil ASA
// 
//  ResInsight is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  ResInsight is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.
// 
//  See the GNU General Public License at <http://www.gnu.org/licenses/gpl.html> 
//  for more details.
//
/////////////////////////////////////////////////////////////////////////////////
#pragma once

#include "cafTensor3.h"
#include "cvfVector3.h"

#include <array>
#include <vector>

class RigFemPart;

class RigFemPartTensorTools
{
public:
    static void calculateElementTensors(const RigFemPart&              part,
                                        const std::vector<caf::Ten3f>& vertexTensors,
                                        std::vector<caf::Ten3f>*       elmTensors);

    static caf::Ten3f calculateElementTensor(const RigFemPart&              part,
                                             const std::vector<caf::Ten3f>& vertexTensors,
                                             size_t                         elementIdx);                                       

    static void calculatePrincipalsAndDirections(const std::vector<caf::Ten3f>&          tensors,
                                                 std::array<std::vector<float>, 3>*      principals,
                                                 std::vector<std::array<cvf::Vec3f, 3>>* principalDirections);
};