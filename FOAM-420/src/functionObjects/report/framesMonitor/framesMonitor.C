/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "framesMonitor/framesMonitor.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "basicThermo/basicThermo.H"
#include "fluidThermo/fluidThermo.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(framesMonitor, 0);
    addToRunTimeSelectionTable(functionObject, framesMonitor, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::framesMonitor::framesMonitor
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    frameName_(dict.lookup<word>("referenceFrame")),
    coorFramePtr_(obr_.lookupObjectRefPtr<coordinateFrame>(frameName_))
{
    framesMonitor::read(dict);
    framesMonitor::writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::framesMonitor::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    return true;
}


bool Foam::functionObjects::framesMonitor::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl << endl;

    // Log to the file
    file() << obr_.time().timeName() << tab
           << coorFramePtr_->transformation() << tab
           << coorFramePtr_->velocity().first() << tab
           << coorFramePtr_->velocity().second() << tab
           << coorFramePtr_->acceleration().first() << tab
           << coorFramePtr_->acceleration().second() << tab
           << endl;

    return true;
}


bool Foam::functionObjects::framesMonitor::write()
{
    return true;
}


void Foam::functionObjects::framesMonitor::writeFileHeader
(
    Ostream& os
) const
{
    Log << "    Logging surface statistics to file: " << fileName_
        << endl;

    writeCommented(os, "Time");
    writeDelimited(os, "transformation");
    writeDelimited(os, "velocity");
    writeDelimited(os, "omega");
    writeDelimited(os, "acceleration");
    writeDelimited(os, "angularAcceleration");

    os << endl;
}

// ************************************************************************* //
