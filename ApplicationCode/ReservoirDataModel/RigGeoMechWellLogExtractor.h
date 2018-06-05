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

#pragma once

#include "RigWellLogExtractor.h"

#include "cafTensor3.h"

#include "cvfBase.h"
#include "cvfObject.h"
#include "cvfMath.h"
#include "cvfStructGrid.h"
#include "cvfVector3.h"

#include <vector>

enum RigElementType;
enum RigFemResultPosEnum;
class RigFemResultAddress;
class RigGeoMechCaseData;
class RigWellPath;

namespace cvf {
    class BoundingBox;
}

//==================================================================================================
/// 
//==================================================================================================
class RigGeoMechWellLogExtractor : public RigWellLogExtractor
{
public:
    RigGeoMechWellLogExtractor(RigGeoMechCaseData* aCase, const RigWellPath* wellpath, const std::string& wellCaseErrorMsgName);

    void                         curveData(const RigFemResultAddress& resAddr, int frameIndex, std::vector<double>* values );
    void                         fractureGradient(const RigFemResultAddress& resAddr, int frameIndex, double rkbDiff, std::vector<double>* values);
    const RigGeoMechCaseData*    caseData()     { return m_caseData.p();}

private:
    class BoreHoleStressCalculator
    {
    public:
        BoreHoleStressCalculator(const caf::Ten3f& tensor, float porePressure, float poissonRatio, float uniaxialCompressiveStrength, int nThetaSubSamples);
        float solveFractureGradient(float minPw, float maxPw, float* thetaOut = nullptr);
        float solveStassiDalia(float minPw, float maxPw, float* thetaOut = nullptr);
        cvf::Vec3f principleStressesAtWall(float pw, float theta) const;
    private:
        typedef float (BoreHoleStressCalculator::*MemberFunc)(float pw, float* thetaOut) const;
        float solveBisection(float minPw, float maxPw, MemberFunc fn, float* thetaOut);
        float solveSecant(MemberFunc fn, float* thetaOut);
        float sigmaTMinOfMin(float wellPressure, float* thetaAtMin) const;
        float stassiDalia(float wellPressure, float* thetaAtMin) const;
        void calculateStressComponents();
        cvf::Vec4f calculateStressComponentsForSegmentAngle(float theta) const;

        caf::Ten3f m_tensor;
        float m_porePressure;
        float m_poissonRatio;
        float m_uniaxialCompressiveStrength;
        int m_nThetaSubSamples;
        std::vector<cvf::Vec4f> m_stressComponents;
    };
   
    template<typename T>
    T                            interpolatedResultValue(RigFemResultPosEnum resultPosType, const std::vector<T>& resultValues, int64_t cpIdx) const;
    void                         calculateIntersection();
    std::vector<size_t>          findCloseCells(const cvf::BoundingBox& bb);
    virtual cvf::Vec3d           calculateLengthInCell(size_t cellIndex, 
                                                       const cvf::Vec3d& startPoint, 
                                                       const cvf::Vec3d& endPoint) const override;
    cvf::Vec3d                   calculateWellPathTangent(int64_t cpIdx) const;
    static caf::Ten3f            transformTensorToWellPathOrientation(const cvf::Vec3d& wellPathTangent,
                                                                      const caf::Ten3f& wellPathTensor);

    cvf::ref<RigGeoMechCaseData> m_caseData;
};


