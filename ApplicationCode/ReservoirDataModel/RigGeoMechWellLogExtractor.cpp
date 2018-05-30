/////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (C) Statoil ASA
//  Copyright (C) Ceetron Solutions AS
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

//==================================================================================================
/// 
//==================================================================================================
#include "RigGeoMechWellLogExtractor.h"
#include "RigFemPart.h"
#include "RigFemPartCollection.h"
#include "RigGeoMechCaseData.h"
#include "RigFemPartResultsCollection.h"

#include "RigWellLogExtractionTools.h"
#include "RigWellPath.h"
#include "RigWellPathIntersectionTools.h"

#include "cafTensor3.h"
#include "cvfGeometryTools.h"


//==================================================================================================
/// Internal root finding class to find a well Pressure that gives a zero SigmaT
//==================================================================================================

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RigGeoMechWellLogExtractor::SigmaCalculator::SigmaCalculator(const caf::Ten3f& tensor, float porePressure, float poissonRatio, int nThetaSubSamples)
    : m_tensor(tensor)
    , m_porePressure(porePressure)
    , m_poissonRatio(poissonRatio)
    , m_nThetaSubSamples(nThetaSubSamples)
{
    calculateStressComponents();
}


//--------------------------------------------------------------------------------------------------
/// Simple bisection method for now
//--------------------------------------------------------------------------------------------------
float RigGeoMechWellLogExtractor::SigmaCalculator::solveForPwBisection(float minPw, float maxPw) const
{
    const int N = 10000;
    const float epsilon = 1.0e-6f;

    float minPwSigma = minimumSigmaT(minPw);
    float maxPwSigma = minimumSigmaT(maxPw);
    float range = maxPw - minPw;
    // Fallback to expand range in case the root is not within the original range
    for (int i = 0; i < 10 && minPwSigma * maxPwSigma > 0.0; ++i)
    {
        minPw -= range;
        maxPw += range;
        range = maxPw - minPw;
        minPwSigma = minimumSigmaT(minPw);
        maxPwSigma = minimumSigmaT(maxPw);
    }

    CVF_ASSERT(minPwSigma * maxPwSigma < 0.0);

    // Bi-section root finding method: https://en.wikipedia.org/wiki/Bisection_method
    int i = 0;
    for (; i <= N && range > maxPw * epsilon; ++i)
    {
        float midPw = (minPw + maxPw) * 0.5;
        float midPwSigma = minimumSigmaT(midPw);
        if (midPwSigma * minPwSigma < 0.0)
        {
            maxPw = midPw;
            maxPwSigma = midPwSigma;
        }
        else
        {
            minPw = midPw;
            minPwSigma = midPwSigma;
        }
        range = maxPw - minPw;
    }
    CVF_ASSERT(i < N); // Otherwise it hasn't converged
    return 0.5 * (minPw + maxPw);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
float RigGeoMechWellLogExtractor::SigmaCalculator::minimumSigmaT(float wellPressure) const
{
    float minSigma = std::numeric_limits<float>::max();
    for (const cvf::Vec3f& stressComponentsForAngle : m_stressComponents)
    {
        float sigma_theta = stressComponentsForAngle[0] - wellPressure;
        const float& sigma_z = stressComponentsForAngle[1];
        float tauSqrx4 = std::pow(stressComponentsForAngle[2], 2) * 4.0;
        float sigma_t = 0.5 * ((sigma_z + sigma_theta) - std::sqrt(std::pow(sigma_z - sigma_theta, 2) + tauSqrx4));
        minSigma = std::min(minSigma, sigma_t);
    }
    return minSigma - m_porePressure;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RigGeoMechWellLogExtractor::SigmaCalculator::calculateStressComponents()
{
    m_stressComponents.reserve(m_nThetaSubSamples);

    for (int i = 0; i < m_nThetaSubSamples; ++i)
    {
        float theta = (i *cvf::PI_F) / (m_nThetaSubSamples - 1.0f);
        cvf::Vec3f stressComponentsForAngle = calculateStressComponentsForSegmentAngle(theta);
        m_stressComponents.push_back(stressComponentsForAngle);
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
cvf::Vec3f RigGeoMechWellLogExtractor::SigmaCalculator::calculateStressComponentsForSegmentAngle(float theta)
{
    cvf::Vec3f stressComponents;

    const float& sx = m_tensor[caf::Ten3f::SXX];
    const float& sy = m_tensor[caf::Ten3f::SYY];
    const float& sz = m_tensor[caf::Ten3f::SZZ];
    const float& txy = m_tensor[caf::Ten3f::SXY];
    const float& txz = m_tensor[caf::Ten3f::SZX];
    const float& tyz = m_tensor[caf::Ten3f::SYZ];

    stressComponents[0] = sx + sy - 2 * (sx - sy) * cos(2 * theta) - 4 * txy * sin(2 * theta);
    stressComponents[1] = sz - m_poissonRatio * (2 * (sx - sy) * cos(2 * theta) + 4 * txy * sin(2 * theta));
    stressComponents[2] = 2 * (tyz * cos(theta) - txz * sin(theta));

    return stressComponents;
}


//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RigGeoMechWellLogExtractor::RigGeoMechWellLogExtractor(RigGeoMechCaseData* aCase,
                                                       const RigWellPath*  wellpath,
                                                       const std::string&  wellCaseErrorMsgName)
    : RigWellLogExtractor(wellpath, wellCaseErrorMsgName)
    , m_caseData(aCase)
{
    calculateIntersection();
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
void RigGeoMechWellLogExtractor::curveData(const RigFemResultAddress& resAddr, int frameIndex, std::vector<double>* values)
{   
    CVF_TIGHT_ASSERT(values);
    
    if (!resAddr.isValid()) return ;

    CVF_ASSERT(resAddr.fieldName != "FracGrad");

    RigFemResultAddress convResAddr = resAddr;

    // When showing POR results, always use the element nodal result, 
    // to get correct handling of elements without POR results
     
    if (convResAddr.fieldName == "POR-Bar") convResAddr.resultPosType = RIG_ELEMENT_NODAL;

    const RigFemPart* femPart                 = m_caseData->femParts()->part(0);
    const std::vector<float>& resultValues    = m_caseData->femPartResults()->resultValues(convResAddr, 0, frameIndex);

    if (!resultValues.size()) return;

    values->resize(m_intersections.size());

    for (size_t cpIdx = 0; cpIdx < m_intersections.size(); ++cpIdx)
    {
        (*values)[cpIdx] = static_cast<double>(interpolatedResultValue<float>(convResAddr, resultValues, cpIdx));
    }
  
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RigGeoMechWellLogExtractor::fractionGradient(const RigFemResultAddress& resAddr, int frameIndex, double rkbDiff, std::vector<double>* values)
{
    CVF_TIGHT_ASSERT(values);
    CVF_ASSERT(resAddr.resultPosType == RIG_INTEGRATION_POINT);
    CVF_ASSERT(resAddr.fieldName == "FracGrad");

    const float poissonRatio = 0.25f; // TODO: Read this in some way

    const RigFemPart* femPart = m_caseData->femParts()->part(0);
    const std::vector<cvf::Vec3f>& nodeCoords = femPart->nodes().coordinates;
    RigFemPartResultsCollection* resultCollection = m_caseData->femPartResults();

    RigFemResultAddress stressResAddr(RIG_INTEGRATION_POINT, std::string("ST"), "");
    stressResAddr.fieldName = std::string("ST");

    RigFemResultAddress porBarResAddr(RIG_ELEMENT_NODAL, std::string("POR-Bar"), "");

    std::vector<caf::Ten3f> vertexStresses = resultCollection->tensors(stressResAddr, 0, frameIndex);
    if (!vertexStresses.size()) return;

    values->resize(m_intersections.size(), 0.0);

    std::vector<float> porePressures = resultCollection->resultValues(porBarResAddr, 0, frameIndex);

//#pragma omp parallel for
    for (int64_t cpIdx = 0; cpIdx < (int64_t) m_intersections.size(); ++cpIdx)
    {        
        size_t elmIdx = m_intersectedCellsGlobIdx[cpIdx];
        RigElementType elmType = femPart->elementType(elmIdx);

        if (!(elmType == HEX8 || elmType == HEX8P)) continue;

        caf::Ten3f interpolatedStress;
        if (cpIdx == 0 || cpIdx == (int64_t)m_intersections.size() - 1)
        {
            interpolatedStress = interpolatedResultValue(stressResAddr, vertexStresses, cpIdx);
        }
        else if (cpIdx % 2 == 0)
        {
            interpolatedStress = (interpolatedResultValue(stressResAddr, vertexStresses, cpIdx - 1)
                               +  interpolatedResultValue(stressResAddr, vertexStresses, cpIdx))  * 0.5;
        }
        else
        {
            interpolatedStress = (interpolatedResultValue(stressResAddr, vertexStresses, cpIdx)
                               + interpolatedResultValue(stressResAddr, vertexStresses, cpIdx + 1))  * 0.5;
        }
        double trueVerticalDepth = -m_intersections[cpIdx].z();
        float porePressure = trueVerticalDepth * 9.81 / 100.0;
        if (!porePressures.empty())
        {
            float interpolatedPorePressure = interpolatedResultValue(porBarResAddr, porePressures, cpIdx);
            if (interpolatedPorePressure != std::numeric_limits<float>::infinity() &&
                interpolatedPorePressure != -std::numeric_limits<float>::infinity())
            {
                porePressure = interpolatedPorePressure;
            }
        }
        
        cvf::Vec3d wellPathTangent = calculateWellPathTangent(cpIdx);
        caf::Ten3f wellPathOrientedStress = transformTensorToWellPathOrientation(wellPathTangent, interpolatedStress);

        SigmaCalculator sigmaCalculator(wellPathOrientedStress, porePressure, poissonRatio, 16);
        float pw_fg = sigmaCalculator.solveForPwBisection(0.0, 2*porePressure);
        float fractureGradient = 100.0 * pw_fg / ((trueVerticalDepth + rkbDiff) * 9.81f);
        (*values)[cpIdx] = fractureGradient;
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
template<typename T>
T RigGeoMechWellLogExtractor::interpolatedResultValue(const RigFemResultAddress& resAddr, const std::vector<T>& resultValues, size_t cpIdx) const
{
    const RigFemPart* femPart = m_caseData->femParts()->part(0);
    const std::vector<cvf::Vec3f>& nodeCoords = femPart->nodes().coordinates;

    size_t elmIdx = m_intersectedCellsGlobIdx[cpIdx];
    RigElementType elmType = femPart->elementType(elmIdx);

    if (!(elmType == HEX8 || elmType == HEX8P)) return T();


    if (resAddr.resultPosType == RIG_ELEMENT)
    {
        return resultValues[elmIdx];        
    }

    cvf::StructGridInterface::FaceType cellFace = m_intersectedCellFaces[cpIdx];

    if (cellFace == cvf::StructGridInterface::NO_FACE)
    {
        // TODO: Should interpolate within the whole hexahedron. This requires hexahedral barycentric calculation.
        // For now just pick the average value for the cell.
        T sumOfVertexValues = resultValues[femPart->elementNodeResultIdx(static_cast<int>(elmIdx), 0)];
        for (int i = 1; i < 8; ++i)
        {
            sumOfVertexValues = sumOfVertexValues + resultValues[femPart->elementNodeResultIdx(static_cast<int>(elmIdx), i)];
        }
        return sumOfVertexValues * (1.0 / 8.0);
    }

    int faceNodeCount = 0;
    const int* faceLocalIndices = RigFemTypes::localElmNodeIndicesForFace(elmType, cellFace, &faceNodeCount);
    const int* elmNodeIndices = femPart->connectivities(elmIdx);

    cvf::Vec3d v0(nodeCoords[elmNodeIndices[faceLocalIndices[0]]]);
    cvf::Vec3d v1(nodeCoords[elmNodeIndices[faceLocalIndices[1]]]);
    cvf::Vec3d v2(nodeCoords[elmNodeIndices[faceLocalIndices[2]]]);
    cvf::Vec3d v3(nodeCoords[elmNodeIndices[faceLocalIndices[3]]]);

    size_t resIdx0 = cvf::UNDEFINED_SIZE_T;
    size_t resIdx1 = cvf::UNDEFINED_SIZE_T;
    size_t resIdx2 = cvf::UNDEFINED_SIZE_T;
    size_t resIdx3 = cvf::UNDEFINED_SIZE_T;

    if (resAddr.resultPosType == RIG_NODAL)
    {
        resIdx0 = elmNodeIndices[faceLocalIndices[0]];
        resIdx1 = elmNodeIndices[faceLocalIndices[1]];
        resIdx2 = elmNodeIndices[faceLocalIndices[2]];
        resIdx3 = elmNodeIndices[faceLocalIndices[3]];
    }
    else
    {
        resIdx0 = (size_t)femPart->elementNodeResultIdx((int)elmIdx, faceLocalIndices[0]);
        resIdx1 = (size_t)femPart->elementNodeResultIdx((int)elmIdx, faceLocalIndices[1]);
        resIdx2 = (size_t)femPart->elementNodeResultIdx((int)elmIdx, faceLocalIndices[2]);
        resIdx3 = (size_t)femPart->elementNodeResultIdx((int)elmIdx, faceLocalIndices[3]);
    }

    T interpolatedValue = cvf::GeometryTools::interpolateQuad<T>(
        v0, resultValues[resIdx0],
        v1, resultValues[resIdx1],
        v2, resultValues[resIdx2],
        v3, resultValues[resIdx3],
        m_intersections[cpIdx]
    );

    return interpolatedValue;
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
void RigGeoMechWellLogExtractor::calculateIntersection()
{
    CVF_ASSERT(m_caseData->femParts()->partCount() == 1);

    std::map<RigMDCellIdxEnterLeaveKey, HexIntersectionInfo > uniqueIntersections;

    const RigFemPart* femPart = m_caseData->femParts()->part(0);
    const std::vector<cvf::Vec3f>& nodeCoords =  femPart->nodes().coordinates;

    for (size_t wpp = 0; wpp < m_wellPath->m_wellPathPoints.size() - 1; ++wpp)
    {
        std::vector<HexIntersectionInfo> intersections;
        cvf::Vec3d p1 = m_wellPath->m_wellPathPoints[wpp];
        cvf::Vec3d p2 = m_wellPath->m_wellPathPoints[wpp+1];

        cvf::BoundingBox bb;

        bb.add(p1);
        bb.add(p2);

        std::vector<size_t> closeCells = findCloseCells(bb);

        cvf::Vec3d hexCorners[8];
        for (size_t ccIdx = 0; ccIdx < closeCells.size(); ++ccIdx)
        {
            RigElementType elmType = femPart->elementType(closeCells[ccIdx]);
            if (!(elmType == HEX8 || elmType == HEX8P)) continue;

            const int* cornerIndices = femPart->connectivities(closeCells[ccIdx]);

            hexCorners[0] = cvf::Vec3d(nodeCoords[cornerIndices[0]]);
            hexCorners[1] = cvf::Vec3d(nodeCoords[cornerIndices[1]]);
            hexCorners[2] = cvf::Vec3d(nodeCoords[cornerIndices[2]]);
            hexCorners[3] = cvf::Vec3d(nodeCoords[cornerIndices[3]]);
            hexCorners[4] = cvf::Vec3d(nodeCoords[cornerIndices[4]]);
            hexCorners[5] = cvf::Vec3d(nodeCoords[cornerIndices[5]]);
            hexCorners[6] = cvf::Vec3d(nodeCoords[cornerIndices[6]]);
            hexCorners[7] = cvf::Vec3d(nodeCoords[cornerIndices[7]]);

            //int intersectionCount = RigHexIntersector::lineHexCellIntersection(p1, p2, hexCorners, closeCells[ccIdx], &intersections);
            RigHexIntersectionTools::lineHexCellIntersection(p1, p2, hexCorners, closeCells[ccIdx], &intersections);
        }

        // Now, with all the intersections of this piece of line, we need to 
        // sort them in order, and set the measured depth and corresponding cell index

        // Inserting the intersections in this map will remove identical intersections
        // and sort them according to MD, CellIdx, Leave/enter

        double md1 = m_wellPath->m_measuredDepths[wpp];
        double md2 = m_wellPath->m_measuredDepths[wpp+1];

        insertIntersectionsInMap(intersections,
                                 p1, md1, p2, md2,
                                 &uniqueIntersections);
    }

    this->populateReturnArrays(uniqueIntersections);
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
std::vector<size_t> RigGeoMechWellLogExtractor::findCloseCells(const cvf::BoundingBox& bb)
{
    std::vector<size_t> closeCells;

    if (m_caseData->femParts()->partCount())
    {
        m_caseData->femParts()->part(0)->findIntersectingCells(bb, &closeCells);
    }
    return closeCells;
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
cvf::Vec3d RigGeoMechWellLogExtractor::calculateLengthInCell(size_t cellIndex, const cvf::Vec3d& startPoint, const cvf::Vec3d& endPoint) const
{
    std::array<cvf::Vec3d, 8> hexCorners;

    const RigFemPart* femPart = m_caseData->femParts()->part(0);
    const std::vector<cvf::Vec3f>& nodeCoords =  femPart->nodes().coordinates;
    const int* cornerIndices = femPart->connectivities(cellIndex);

    hexCorners[0] = cvf::Vec3d(nodeCoords[cornerIndices[0]]);
    hexCorners[1] = cvf::Vec3d(nodeCoords[cornerIndices[1]]);
    hexCorners[2] = cvf::Vec3d(nodeCoords[cornerIndices[2]]);
    hexCorners[3] = cvf::Vec3d(nodeCoords[cornerIndices[3]]);
    hexCorners[4] = cvf::Vec3d(nodeCoords[cornerIndices[4]]);
    hexCorners[5] = cvf::Vec3d(nodeCoords[cornerIndices[5]]);
    hexCorners[6] = cvf::Vec3d(nodeCoords[cornerIndices[6]]);
    hexCorners[7] = cvf::Vec3d(nodeCoords[cornerIndices[7]]);

    return RigWellPathIntersectionTools::calculateLengthInCell(hexCorners, startPoint, endPoint); 
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
cvf::Vec3d RigGeoMechWellLogExtractor::calculateWellPathTangent(int64_t cpIdx) const
{
    cvf::Vec3d wellPathTangent;
    if (cpIdx % 2 == 0)
    {
        wellPathTangent = m_intersections[cpIdx + 1] - m_intersections[cpIdx];
    }
    else
    {
        wellPathTangent = m_intersections[cpIdx] - m_intersections[cpIdx - 1];
    }
    CVF_ASSERT(wellPathTangent.length() > 1.0e-7);
    return wellPathTangent.getNormalized();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
caf::Ten3f RigGeoMechWellLogExtractor::transformTensorToWellPathOrientation(const cvf::Vec3d& wellPathTangent,
                                                                            const caf::Ten3f& tensor)
{
    // Create local coordinate system for well path segment
    cvf::Vec3f local_z = cvf::Vec3f(wellPathTangent);
    cvf::Vec3f local_x = local_z.perpendicularVector().getNormalized();
    cvf::Vec3f local_y = (local_z ^ local_x).getNormalized();
    // Calculate the rotation matrix from global i, j, k to local x, y, z.
    cvf::Mat4f rotationMatrix = cvf::Mat4f::fromCoordSystemAxes(&local_x, &local_y, &local_z);

    return tensor.rotated(rotationMatrix.toMatrix3());
}

