#include "RigFemPartTensorTools.h"

#include "RigFemPart.h"


//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RigFemPartTensorTools::calculateElementTensors(const RigFemPart&              part,
                                                    const std::vector<caf::Ten3f>& vertexTensors,
                                                    std::vector<caf::Ten3f>*       elmTensors)
{
    CVF_ASSERT(elmTensors);

    size_t elmCount = part.elementCount();
    elmTensors->resize(elmCount);

    for (size_t elmIdx = 0; elmIdx < elmCount; elmIdx++)
    {
        if (RigFemTypes::elmentNodeCount(part.elementType(elmIdx)) == 8)
        {
            (*elmTensors)[elmIdx] = calculateElementTensor(part, vertexTensors, elmIdx);
        }
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
caf::Ten3f RigFemPartTensorTools::calculateElementTensor(const RigFemPart& part, const std::vector<caf::Ten3f>& vertexTensors, size_t elmIdx)
{
    CVF_ASSERT(RigFemTypes::elmentNodeCount(part.elementType(elmIdx)) == 8);    
    caf::Ten3f tensorSumOfElmNodes = vertexTensors[part.elementNodeResultIdx(static_cast<int>(elmIdx), 0)];
    for (int i = 1; i < 8; i++)
    {
        tensorSumOfElmNodes = tensorSumOfElmNodes + vertexTensors[part.elementNodeResultIdx(static_cast<int>(elmIdx), i)];
    }
    return tensorSumOfElmNodes * (1.0 / 8.0);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RigFemPartTensorTools::calculatePrincipalsAndDirections(const std::vector<caf::Ten3f>&          tensors,
                                                             std::array<std::vector<float>, 3>*      principals,
                                                             std::vector<std::array<cvf::Vec3f, 3>>* principalDirections)
{
    CVF_ASSERT(principals);
    CVF_ASSERT(principalDirections);

    size_t elmCount = tensors.size();

    (*principals)[0].resize(elmCount);
    (*principals)[1].resize(elmCount);
    (*principals)[2].resize(elmCount);

    (*principalDirections).resize(elmCount);

    for (size_t nIdx = 0; nIdx < elmCount; ++nIdx)
    {
        cvf::Vec3f principalDirs[3];
        cvf::Vec3f principalValues = tensors[nIdx].calculatePrincipals(principalDirs);

        (*principals)[0][nIdx] = principalValues[0];
        (*principals)[1][nIdx] = principalValues[1];
        (*principals)[2][nIdx] = principalValues[2];

        (*principalDirections)[nIdx][0] = principalDirs[0];
        (*principalDirections)[nIdx][1] = principalDirs[1];
        (*principalDirections)[nIdx][2] = principalDirs[2];
    }
}