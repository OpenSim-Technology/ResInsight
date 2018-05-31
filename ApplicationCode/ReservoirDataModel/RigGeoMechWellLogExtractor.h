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
    class SigmaCalculator
    {
    public:
        SigmaCalculator(const caf::Ten3f& tensor, float porePressure, float poissonRatio, int nThetaSubSamples);
        void solveForPwBisection(float minPw, float maxPw);
        cvf::Vec3f principalStressesPlusPorePressure() const;
    private:
        std::pair<float, float> sigmaTMinMax() const;
        float sigmaTMinOfMin(float wellPressure, float* thetaAtMin) const;
        void calculateStressComponents();
        cvf::Vec4f calculateStressComponentsForSegmentAngle(float theta) const;

        caf::Ten3f m_tensor;
        float m_porePressure;
        float m_poissonRatio;
        int m_nThetaSubSamples;
        std::vector<cvf::Vec4f> m_stressComponents;

        float m_wellPressureFractureGradient;
        float m_thetaPrincipleStress;
    };

    cvf::Vec3f calculatePrincipalStressesPlusPorePressure(RigFemResultPosEnum            resultPosType,
                                                          const std::vector<caf::Ten3f>& vertexStresses,
                                                          const std::vector<float>&      porePressures,
                                                          int64_t                        cpIdx) const;
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


