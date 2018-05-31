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
RigGeoMechWellLogExtractor::BoreHoleStressCalculator::BoreHoleStressCalculator(const caf::Ten3f& tensor,
                                                             float             porePressure,
                                                             float             poissonRatio,
                                                             float             uniaxialCompressiveStrength,
                                                             int               nThetaSubSamples)
    : m_tensor(tensor)
    , m_porePressure(porePressure)
    , m_poissonRatio(poissonRatio)
    , m_uniaxialCompressiveStrength(uniaxialCompressiveStrength)
    , m_nThetaSubSamples(nThetaSubSamples)
{
    calculateStressComponents();
}


//--------------------------------------------------------------------------------------------------
/// Simple bisection method for now
//--------------------------------------------------------------------------------------------------
float RigGeoMechWellLogExtractor::BoreHoleStressCalculator::solveFractureGradient(float minPw, float maxPw, float* thetaOut)
{
    MemberFunc fn = &RigGeoMechWellLogExtractor::BoreHoleStressCalculator::sigmaTMinOfMin;
    return solveSecant(fn, thetaOut);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
float RigGeoMechWellLogExtractor::BoreHoleStressCalculator::solveStassiDalia(float minPw, float maxPw, float* thetaOut)
{
    MemberFunc fn = &RigGeoMechWellLogExtractor::BoreHoleStressCalculator::stassiDalia;
    return solveSecant(fn, thetaOut);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
cvf::Vec3f RigGeoMechWellLogExtractor::BoreHoleStressCalculator::principleStressesAtWall(float pw, float theta) const
{
    cvf::Vec4f stressComponentsForAngle = calculateStressComponentsForSegmentAngle(theta);
    float sigma_theta = stressComponentsForAngle[1] - pw;
    const float& sigma_z = stressComponentsForAngle[2];
    float tauSqrx4 = std::pow(stressComponentsForAngle[3], 2) * 4.0;

    float sigmaComponent1 = sigma_z + sigma_theta;
    float sigmaComponent2 = std::sqrt(std::pow(sigma_z - sigma_theta, 2) + tauSqrx4);    
    return cvf::Vec3f(pw - m_porePressure, 0.5 * (sigmaComponent1 + sigmaComponent2) - m_porePressure, 0.5 * (sigmaComponent1 - sigmaComponent2) - m_porePressure);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
float RigGeoMechWellLogExtractor::BoreHoleStressCalculator::solveBisection(float minPw, float maxPw, MemberFunc fn, float* thetaOut)
{
    const int N = 50;
    const float epsilon = 1.0e-6f;

    float theta = 0.0;

    float minPwFuncVal = std::invoke(fn, this, minPw, &theta);
    float maxPwFuncVal = std::invoke(fn, this, maxPw, &theta);
    float range = maxPw - minPw;
    
    // Bi-section root finding method: https://en.wikipedia.org/wiki/Bisection_method
    int i = 0;
    for (; i <= N && range > m_porePressure * epsilon; ++i)
    {
        float midPw = (minPw + maxPw) * 0.5;
        float midPwFuncVal = std::invoke(fn, this, midPw, &theta);
        if (midPwFuncVal * minPwFuncVal < 0.0)
        {
            maxPw = midPw;
            maxPwFuncVal = midPwFuncVal;
        }
        else
        {
            minPw = midPw;
            minPwFuncVal = midPwFuncVal;
        }
        range = maxPw - minPw;
    }
    CVF_ASSERT(i < N); // Otherwise it hasn't converged
                       // Return linear solution between minPw and maxPw.

    if (thetaOut)
    {
        *thetaOut = theta;
    }

    return minPw - minPwFuncVal * (maxPw - minPw) / (maxPwFuncVal - minPwFuncVal);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
float RigGeoMechWellLogExtractor::BoreHoleStressCalculator::solveRegularFalsi(MemberFunc fn, float* thetaOut)
{
    const int N = 100;
    const float epsilon = 1.0e-6f;

    float theta = 0.0;

    float xl = 0.0;
    float xr = m_porePressure;

    float f_xl = std::invoke(fn, this, xl, &theta);
    float f_xr = std::invoke(fn, this, xr, &theta);

    float xm = xm = (xl * f_xr - xr * f_xl) / (f_xr - f_xl);
    float f_xm = std::invoke(fn, this, xm, &theta);

    if (f_xr * f_xl > 0.0) // Original interval was bad
    {
        if (f_xm * f_xl < 0.0)
        {
            xr = xm;
            f_xr = f_xm;
        }
        else if (f_xm * f_xr < 0.0)
        {
            xl = xm;
            f_xl = f_xm;
        }
        else
        {
            CVF_ASSERT(false); // Bad interval
        }
    }


    // Regular Falsi root finding method: https://en.wikipedia.org/wiki/False_position_method
    int i = 0;
    int side = 0;
    for (; i <= N; ++i)
    {
        xm = (xl * f_xr - xr * f_xl) / (f_xr - f_xl);
        if (std::abs(xr - xl) < epsilon * std::abs(xr + xl)) break;
        float f_xm = std::invoke(fn, this, xm, &theta);

        if (f_xm * f_xr > 0.0)
        {
            xr = xm;
            f_xr = f_xm;
            if (side == -1) f_xl /= 2;
            side = -1;
        }
        else if (f_xm * f_xl > 0.0)
        {
            xl = xm;
            f_xl = f_xm;
            if (side == 1) f_xr /= 2;
            side = 1;
        }
        else {
            break;
        }
    }
    CVF_ASSERT(i < N); // Otherwise it hasn't converged
                       // Return linear solution between minPw and maxPw.

    if (thetaOut)
    {
        *thetaOut = theta;
    }

    return xr;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
float RigGeoMechWellLogExtractor::BoreHoleStressCalculator::solveSecant(MemberFunc fn, float* thetaOut)
{
    const float epsilon = 1.0e-8f;
    const int N = 200;
    float theta = 0.0;

    float x_0 = 0.0;    
    float f_x0 = std::invoke(fn, this, x_0, &theta);
    float x_1 = m_porePressure;
    float f_x1 = std::invoke(fn, this, x_1, &theta);
    float x = 0.0;
    int i = 0;
    for (; i < N && std::abs(f_x1 - f_x0) > m_porePressure * epsilon; ++i)
    {
        x = x_1 - f_x1 * (x_1 - x_0) / (f_x1 - f_x0);
        float f_x = std::invoke(fn, this, x, &theta);

        // Update iteration variables
        x_0 = x_1;
        f_x0 = f_x1;
        x_1 = x;
        f_x1 = f_x;        
    }
    CVF_ASSERT(i < N); // Otherwise it hasn't converged
                       // Return linear solution between minPw and maxPw.

    if (thetaOut)
    {
        *thetaOut = theta;
    }
    return x;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
float RigGeoMechWellLogExtractor::BoreHoleStressCalculator::sigmaTMinOfMin(float wellPressure, float* thetaAtMin) const
{
    CVF_ASSERT(thetaAtMin);
    float sigma_t_min_min = std::numeric_limits<float>::max();
    for (const cvf::Vec4f& stressComponentsForAngle : m_stressComponents)
    {
        float sigma_theta = stressComponentsForAngle[1] - wellPressure;
        const float& sigma_z = stressComponentsForAngle[2];
        float tauSqrx4 = std::pow(stressComponentsForAngle[3], 2) * 4.0;
        float sigma_t_min = 0.5 * ((sigma_z + sigma_theta) - std::sqrt(std::pow(sigma_z - sigma_theta, 2) + tauSqrx4));
        if (sigma_t_min < sigma_t_min_min)
        {
            sigma_t_min_min = sigma_t_min;
            *thetaAtMin = stressComponentsForAngle[0];
        }
    }
    return sigma_t_min_min - m_porePressure;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
float RigGeoMechWellLogExtractor::BoreHoleStressCalculator::stassiDalia(float wellPressure, float* thetaAtMin) const
{
    CVF_ASSERT(thetaAtMin);
    float minStassiDalia = std::numeric_limits<float>::max();
    for (const cvf::Vec4f& stressComponentsForAngle : m_stressComponents)
    {
        float sigma_theta = stressComponentsForAngle[1] - wellPressure;
        const float& sigma_z = stressComponentsForAngle[2];
        float tauSqrx4 = std::pow(stressComponentsForAngle[3], 2) * 4.0;

        float sigma_1 = wellPressure - m_porePressure;
        float sigma_2 = 0.5 * ((sigma_z + sigma_theta) + std::sqrt(std::pow(sigma_z - sigma_theta, 2) + tauSqrx4)) - m_porePressure;
        float sigma_3 = 0.5 * ((sigma_z + sigma_theta) - std::sqrt(std::pow(sigma_z - sigma_theta, 2) + tauSqrx4)) - m_porePressure;

        float stassiDalia = std::pow(sigma_1 - sigma_2, 2) + std::pow(sigma_2 - sigma_3, 2) + std::pow(sigma_1 - sigma_3, 2)
                          - 2 * m_uniaxialCompressiveStrength * (sigma_1 + sigma_2 + sigma_3);

        if (stassiDalia < minStassiDalia)
        {
            minStassiDalia = stassiDalia;
            *thetaAtMin = stressComponentsForAngle[0];
        }
    }
    return minStassiDalia;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RigGeoMechWellLogExtractor::BoreHoleStressCalculator::calculateStressComponents()
{
    m_stressComponents.reserve(m_nThetaSubSamples);

    for (int i = 0; i < m_nThetaSubSamples; ++i)
    {
        float theta = (i *cvf::PI_F) / (m_nThetaSubSamples - 1.0f);
        cvf::Vec4f stressComponentsForAngle = calculateStressComponentsForSegmentAngle(theta);
        m_stressComponents.push_back(stressComponentsForAngle);
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
cvf::Vec4f RigGeoMechWellLogExtractor::BoreHoleStressCalculator::calculateStressComponentsForSegmentAngle(float theta) const
{
    cvf::Vec4f stressComponents;

    const float& sx = m_tensor[caf::Ten3f::SXX];
    const float& sy = m_tensor[caf::Ten3f::SYY];
    const float& sz = m_tensor[caf::Ten3f::SZZ];
    const float& txy = m_tensor[caf::Ten3f::SXY];
    const float& txz = m_tensor[caf::Ten3f::SZX];
    const float& tyz = m_tensor[caf::Ten3f::SYZ];

    stressComponents[0] = theta;
    stressComponents[1] = sx + sy - 2 * (sx - sy) * cos(2 * theta) - 4 * txy * sin(2 * theta);
    stressComponents[2] = sz - m_poissonRatio * (2 * (sx - sy) * cos(2 * theta) + 4 * txy * sin(2 * theta));
    stressComponents[3] = 2 * (tyz * cos(theta) - txz * sin(theta));

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
        (*values)[cpIdx] = static_cast<double>(interpolatedResultValue<float>(convResAddr.resultPosType, resultValues, cpIdx));
    }
  
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RigGeoMechWellLogExtractor::fractureGradient(const RigFemResultAddress& resAddr, int frameIndex, double rkbDiff, std::vector<double>* values)
{
    const float poissonRatio = 0.25f; // TODO: Read this in
    const float uniaxialStrengthInBars = 20.0f; // TODO: Read this in

    CVF_TIGHT_ASSERT(values);
    CVF_ASSERT(resAddr.fieldName == "FracGrad" || resAddr.fieldName == "StassidAlia");

    bool fracGrad = resAddr.fieldName == "FracGrad";
    const RigFemPart* femPart = m_caseData->femParts()->part(0);
    const std::vector<cvf::Vec3f>& nodeCoords = femPart->nodes().coordinates;
    RigFemPartResultsCollection* resultCollection = m_caseData->femPartResults();

    RigFemResultAddress stressResAddr(resAddr.resultPosType, std::string("ST"), "");
    stressResAddr.fieldName = std::string("ST");

    RigFemResultAddress porBarResAddr(RIG_ELEMENT_NODAL, std::string("POR-Bar"), "");

    std::vector<caf::Ten3f> vertexStresses = resultCollection->tensors(stressResAddr, 0, frameIndex);
    if (!vertexStresses.size()) return;

    values->resize(m_intersections.size(), 0.0f);

    std::vector<float> porePressures = resultCollection->resultValues(porBarResAddr, 0, frameIndex);

//#pragma omp parallel for
    for (int64_t cpIdx = 0; cpIdx < (int64_t) m_intersections.size(); ++cpIdx)
    {        
        size_t elmIdx = m_intersectedCellsGlobIdx[cpIdx];
        RigElementType elmType = femPart->elementType(elmIdx);

        if (!(elmType == HEX8 || elmType == HEX8P)) continue;


        double trueVerticalDepth = -m_intersections[cpIdx].z();
        float porePressure = trueVerticalDepth * 9.81 / 100.0;
        if (!porePressures.empty())
        {
            float interpolatedPorePressure = interpolatedResultValue(porBarResAddr.resultPosType, porePressures, cpIdx);
            if (interpolatedPorePressure != std::numeric_limits<float>::infinity() &&
                interpolatedPorePressure != -std::numeric_limits<float>::infinity())
            {
                porePressure = interpolatedPorePressure;
            }
        }

        caf::Ten3f interpolatedStress = interpolatedResultValue(stressResAddr.resultPosType, vertexStresses, cpIdx);

        cvf::Vec3d wellPathTangent = calculateWellPathTangent(cpIdx);
        caf::Ten3f wellPathOrientedStress = transformTensorToWellPathOrientation(wellPathTangent, interpolatedStress);

        BoreHoleStressCalculator sigmaCalculator(wellPathOrientedStress, porePressure, poissonRatio, uniaxialStrengthInBars, 16);
        float wellPressure = 0.0;
        if (fracGrad)
        {
            wellPressure = sigmaCalculator.solveFractureGradient(0.0, 3 * porePressure);
        }
        else
        {
            wellPressure = sigmaCalculator.solveStassiDalia(0.0, 3 * porePressure);
        }

        (*values)[cpIdx] = wellPressure / ((-m_intersections[cpIdx].z() + rkbDiff) * 9.81 / 100.0f);
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
template<typename T>
T RigGeoMechWellLogExtractor::interpolatedResultValue(RigFemResultPosEnum resultPosType, const std::vector<T>& resultValues, int64_t cpIdx) const
{
    const RigFemPart* femPart = m_caseData->femParts()->part(0);
    const std::vector<cvf::Vec3f>& nodeCoords = femPart->nodes().coordinates;

    size_t elmIdx = m_intersectedCellsGlobIdx[cpIdx];
    RigElementType elmType = femPart->elementType(elmIdx);

    if (!(elmType == HEX8 || elmType == HEX8P)) return T();


    if (resultPosType == RIG_ELEMENT)
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

    std::vector<size_t> nodeResIdx(4, cvf::UNDEFINED_SIZE_T);

    if (resultPosType == RIG_NODAL)
    {
        for (size_t i = 0; i < nodeResIdx.size(); ++i)
        {
            nodeResIdx[i] = elmNodeIndices[faceLocalIndices[i]];
        }
    }
    else
    {
        for (size_t i = 0; i < nodeResIdx.size(); ++i)
        {
            nodeResIdx[i] = (size_t)femPart->elementNodeResultIdx((int)elmIdx, faceLocalIndices[i]);
        }
    }

    std::vector<T> nodeResultValues;
    nodeResultValues.reserve(4);
    if (resultPosType == RIG_ELEMENT_NODAL || resultPosType == RIG_INTEGRATION_POINT)
    {
        for (size_t i = 0; i < nodeResIdx.size(); ++i)
        {
            int nodeIndex = femPart->nodeIdxFromElementNodeResultIdx(nodeResIdx[i]);
            const std::vector<std::pair<int, int>>& elements = femPart->elementsUsingNode(nodeIndex);
            size_t otherNodeResIdx = femPart->elementNodeResultIdx(elements[0].first, elements[0].second);
            T nodeResultValue = resultValues[otherNodeResIdx];
            for (size_t j = 1; j < elements.size(); ++j)
            {
                otherNodeResIdx = femPart->elementNodeResultIdx(elements[j].first, elements[j].second);
                nodeResultValue = nodeResultValue + resultValues[otherNodeResIdx];
            }
            nodeResultValue = nodeResultValue * (1.0 / elements.size());
            nodeResultValues.push_back(nodeResultValue);
        }
    }
    else {
        for (size_t i = 0; i < nodeResIdx.size(); ++i)
        {
            nodeResultValues.push_back(resultValues[nodeResIdx[i]]);
        }
    }

    T interpolatedValue = cvf::GeometryTools::interpolateQuad<T>(
        v0, nodeResultValues[0],
        v1, nodeResultValues[1],
        v2, nodeResultValues[2],
        v3, nodeResultValues[3],
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

