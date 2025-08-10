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
    (c) 2018 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceAreaReport/surfaceAreaReport.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(surfaceAreaReport, 0);
    addToRunTimeSelectionTable(functionObject, surfaceAreaReport, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::surfaceAreaReport::writeFileHeader(Ostream& os) const
{
    Log << "    Logging surface statistics to file: " << fileName_
        << endl;

    writeCommented(os, "Time");

    writeDelimited(os,"surfaceArea");

    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

functionObjects::surfaceAreaReport::surfaceAreaReport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    patchName_("none"),
    patchID_(-1),
    surfaceArea_(0)
{
    read(dict);
    writeFileHeader(file());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool functionObjects::surfaceAreaReport::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    Log << type() << " " << name() <<  " read:" << nl;

    //determine compressible/incompressible
    {
        if (dict.found("patchName"))
        {
            dict.lookup("patchName") >> patchName_;
        }
        else
        {
            WarningInFunction
                << patchName_ << " could not be found in the database, "
                << "deactivating."
                << endl;
        }

    }

    return true;
}


bool Foam::functionObjects::surfaceAreaReport::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    //calculate
    calculate();

    //Write to screen
    Log << "Surface area [m^2] = ";
    Log << surfaceArea_ << endl;

    file() << obr_.time().timeName() << " ";
    file() << surfaceArea_ << " ";

    file() << endl;

    Log << endl;

    return true;
}

void functionObjects::surfaceAreaReport::calculate()
{
    patchID_ = mesh_.boundaryMesh().findPatchID(patchName_);

    surfaceArea_ = 0;

    if (patchID_ == -1)
    {
        WarningInFunction
            << patchName_ << " patch name not existed"
            << endl;
    }
    else
    {
        surfaceArea_ = gSum(mesh_.boundary()[patchID_].magSf());
    }
}

bool functionObjects::surfaceAreaReport::write()
{
    return true;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
