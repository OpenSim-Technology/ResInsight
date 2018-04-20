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

#include "RicPaste3dWellLogExtractionCurve.h"

#include "RicPasteFeatureImpl.h"
#include "Rim3dWellLogCurveCollection.h"
#include "Rim3dWellLogExtractionCurve.h"
#include "RiuMainWindow.h"

#include "cafPdmObjectGroup.h"
#include "cafSelectionManager.h"

#include <QAction>

CAF_CMD_SOURCE_INIT(RicPaste3dWellLogExtractionCurve, "RicPaste3dWellLogExtractionCurve");

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
bool RicPaste3dWellLogExtractionCurve::isCommandEnabled()
{
    caf::PdmObjectGroup objectGroup;
    RicPasteFeatureImpl::findObjectsFromClipboardRefs(&objectGroup);

    std::vector<caf::PdmPointer<Rim3dWellLogExtractionCurve>> wellLogCurveObjects;
    objectGroup.objectsByType(&wellLogCurveObjects);

    if (wellLogCurveObjects.empty() && wellLogCurveBoxObjects.empty())
    {
        return false;
    }

    caf::PdmObjectHandle* destinationObject =
        dynamic_cast<caf::PdmObjectHandle*>(caf::SelectionManager::instance()->selectedItem());

    if (find3dWellLogCurveCollection(destinationObject))
    {
        return true;
    }

    return false;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicPaste3dWellLogExtractionCurve::onActionTriggered(bool isChecked)
{
    caf::PdmObjectHandle* destinationObject =
        dynamic_cast<caf::PdmObjectHandle*>(caf::SelectionManager::instance()->selectedItem());

    Rim3dWellLogCurveCollection* wellLogCurveCollection =
        RicPaste3dWellLogExtractionCurve::findIntersectionCollection(destinationObject);

    CAF_ASSERT(wellLogCurveCollection);

    caf::PdmObjectGroup objectGroup;
    RicPasteFeatureImpl::findObjectsFromClipboardRefs(&objectGroup);

    if (objectGroup.objects.size() == 0) return;

    std::vector<caf::PdmPointer<Rim3dWellLogExtractionCurve>> wellLogCurveObjects;
    objectGroup.objectsByType(&wellLogCurveObjects);

    for (size_t i = 0; i < wellLogCurveObjects.size(); i++)
    {
        Rim3dWellLogExtractionCurve* wellLogCurve = dynamic_cast<Rim3dWellLogExtractionCurve*>(
            wellLogCurveObjects[i]->xmlCapability()->copyByXmlSerialization(caf::PdmDefaultObjectFactory::instance()));

        QString nameOfCopy = QString("Copy of ") + wellLogCurve->name;
        wellLogCurve->name = nameOfCopy;

        if (i == wellLogCurveObjects.size() - 1)
        {
            wellLogCurveCollection->add3dWellLogCurve(wellLogCurve);
            // Last item to paste. Update!
            wellLogCurveCollection->redrawAffectedViewsAndEditor();
        }
        else
        {
            wellLogCurveCollection->add3dWellLogCurve(wellLogCurve);
        }
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicPaste3dWellLogExtractionCurve::setupActionLook(QAction* actionToSetup)
{
    actionToSetup->setText("Paste (3d Well Log Extraction Curve(s))");

    RicPasteFeatureImpl::setIconAndShortcuts(actionToSetup);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
Rim3dWellLogCurveCollection* RicPaste3dWellLogExtractionCurve::find3dWellLogCurveCollection(caf::PdmObjectHandle* objectHandle)
{
    Rim3dWellLogCurveCollection* wellLogCurveCollection = dynamic_cast<Rim3dWellLogCurveCollection*>(objectHandle);
    if (wellLogCurveCollection)
    {
        return wellLogCurveCollection;
    }

    Rim3dWellLogExtractionCurve* wellLogCurve = dynamic_cast<Rim3dWellLogExtractionCurve*>(objectHandle);
    if (wellLogCurve)
    {
        wellLogCurve->firstAncestorOrThisOfType(wellLogCurveCollection);
        return wellLogCurveCollection;
    }

    return nullptr;
}
