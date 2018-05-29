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

#include "Rim3dWellLogCurve.h"

#include "cafPdmChildField.h"
#include "cafPdmField.h"
#include "cafPdmPtrField.h"

class Rim3dView;
class RimCase;
class RimGeoMechResultDefinition;
class RimEclipseResultDefinition;
class RimWellLogExtractionCurveNameConfig;

//==================================================================================================
///
///
//==================================================================================================
class Rim3dWellLogExtractionCurve : public Rim3dWellLogCurve
{
    CAF_PDM_HEADER_INIT;

public:
    Rim3dWellLogExtractionCurve();
    virtual ~Rim3dWellLogExtractionCurve();

    void            setPropertiesFromView(Rim3dView* view);
    virtual QString resultPropertyString() const override;

    virtual bool    useCurrentTimeStep() const override;
    virtual void    curveValuesAndMds(std::vector<double>* values, std::vector<double>* measuredDepthValues) const override;
    virtual void    curveValuesAndMdsAtTimeStep(std::vector<double>* values, std::vector<double>* measuredDepthValues, int timeStep = -1) const override;
    virtual std::pair<double, double> findCurveValueRange() override;

    virtual QString name() const override;
    virtual QString createCurveAutoName() const override;
    double          rkbDiff() const;
protected:
    virtual caf::PdmFieldHandle*            userDescriptionField() override;
    virtual void                            fieldChangedByUi(const caf::PdmFieldHandle* changedField, const QVariant& oldValue, const QVariant& newValue) override;
private:
    virtual QList<caf::PdmOptionItemInfo> calculateValueOptions(const caf::PdmFieldHandle* fieldNeedingOptions,
                                                                bool*                      useOptionsOnly) override;

    virtual void defineUiOrdering(QString uiConfigName, caf::PdmUiOrdering& uiOrdering) override;
    virtual void initAfterRead() override;
    QString      wellDate() const;
private:
    caf::PdmPtrField<RimCase*>                      m_case;
    caf::PdmField<int>                              m_timeStep;

    caf::PdmChildField<RimEclipseResultDefinition*> m_eclipseResultDefinition;
    caf::PdmChildField<RimGeoMechResultDefinition*> m_geomResultDefinition;
    caf::PdmChildField<RimWellLogExtractionCurveNameConfig*> m_nameConfig;
};
