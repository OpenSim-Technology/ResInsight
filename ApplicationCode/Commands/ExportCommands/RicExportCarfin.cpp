/////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (C) 2017 Statoil ASA
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

#include "RicExportCarfin.h"

#include "RiaLogging.h"

#include "RicCellRangeUi.h"
#include "RicExportCarfinUi.h"
#include "RicExportFeatureImpl.h"

#include "RifEclipseDataTableFormatter.h"

#include "RimEclipseCase.h"

#include "cafPdmUiPropertyViewDialog.h"
#include "cafSelectionManager.h"

#include <QAction>
#include <QFile>

CAF_CMD_SOURCE_INIT(RicExportCarfin, "RicExportCarfin");

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
bool RicExportCarfin::isCommandEnabled()
{
    if (RicExportCarfin::selectedCase() != nullptr)
    {
        return true;
    }

    return false;
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
void RicExportCarfin::onActionTriggered(bool isChecked)
{
    RimEclipseCase* rimCase = RicExportCarfin::selectedCase();
    CVF_ASSERT(rimCase);

    RicExportCarfinUi carfinUi;
    carfinUi.setCase(rimCase);

    caf::PdmUiPropertyViewDialog propertyDialog(nullptr, &carfinUi, "Export CARFIN to Eclipse Data", "");
    RicExportFeatureImpl::configureForExport(&propertyDialog);

    if (propertyDialog.exec() == QDialog::Accepted)
    {
        QString filePath = carfinUi.exportFileName();
        QFile exportFile(filePath);
        if (!exportFile.open(QIODevice::WriteOnly))
        {
            RiaLogging::error(QString("Export CARFIN: Could not open the file: %1").arg(filePath));
            return;
        }

        QTextStream stream(&exportFile);
        RifEclipseDataTableFormatter formatter(stream);

        std::vector<RifEclipseOutputTableColumn> header = {
            RifEclipseOutputTableColumn("I1"),
            RifEclipseOutputTableColumn("I2"),
            RifEclipseOutputTableColumn("J1"),
            RifEclipseOutputTableColumn("J2"),
            RifEclipseOutputTableColumn("K1"),
            RifEclipseOutputTableColumn("K2"),
            RifEclipseOutputTableColumn("NX"),
            RifEclipseOutputTableColumn("NY"),
            RifEclipseOutputTableColumn("NZ"),
            RifEclipseOutputTableColumn("NWMAX"),
            RifEclipseOutputTableColumn("Parent LGR")
        };

        formatter.keyword("CARFIN");
        formatter.header(header);

        formatter.add(carfinUi.cellRange()->start().i());
        formatter.add(carfinUi.cellRange()->start().i() + carfinUi.cellRange()->count().i());

        formatter.add(carfinUi.cellRange()->start().j());
        formatter.add(carfinUi.cellRange()->start().j() + carfinUi.cellRange()->count().j());

        formatter.add(carfinUi.cellRange()->start().k());
        formatter.add(carfinUi.cellRange()->start().k() + carfinUi.cellRange()->count().k());

        formatter.add(carfinUi.lgrCellCount().i());
        formatter.add(carfinUi.lgrCellCount().j());
        formatter.add(carfinUi.lgrCellCount().k());

        formatter.add(carfinUi.maxWellCount());

        if (!carfinUi.gridName().isEmpty())
        {
            formatter.add(carfinUi.gridName());
        }

        formatter.rowCompleted();
        formatter.tableCompleted();
    }
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
void RicExportCarfin::setupActionLook(QAction* actionToSetup)
{
    actionToSetup->setText("Export CARFIN ...");
    actionToSetup->setIcon(QIcon(":/Save.png"));
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
RimEclipseCase* RicExportCarfin::selectedCase()
{
    std::vector<RimEclipseCase*> selectedObjects;
    caf::SelectionManager::instance()->objectsByType(&selectedObjects);

    if (selectedObjects.size() == 1)
    {
        return selectedObjects[0];
    }

    return nullptr;
}