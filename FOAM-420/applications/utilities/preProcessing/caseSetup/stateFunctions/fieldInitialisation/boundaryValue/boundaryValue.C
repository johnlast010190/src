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
    (c) 2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "boundaryValue.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(boundaryValue, 0);
    addToRunTimeSelectionTable(fieldInit, boundaryValue, initMethod);
}
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldInitializations::boundaryValue::boundaryValue
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDict,
    const word& fN
)
:
    fieldInit(mesh, localDb, fieldDict, fN)
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fieldInitializations::boundaryValue::correct()
{
    //if already initialised stop
    if (initialised())
    {
        return;
    }

    fieldInit::initMsg();

    word patchName = word
    (
        initDict().lookup("patch")
    );

    label patchID = mesh().boundaryMesh().findPatchID(patchName);

    if (patchID < 0)
    {
        FatalError << "Setup error: patch name " << patchName
                   << " specified as part of `boundaryValue`"
                   << " initialisation of field " << name()
                   << " can not be found." << exit(FatalError);
    }

    Info<< tab << "Using the average value from patch " << patchName << endl;

    if (setBoundaryValue<scalar>(patchID));
    else if (setBoundaryValue<vector>(patchID));
    else if (setBoundaryValue<tensor>(patchID));
    else if (setBoundaryValue<symmTensor>(patchID));
    else if (setBoundaryValue<sphericalTensor>(patchID));
    else
    {
        WarningInFunction
            << "Could not find field " << name()
            << " in database. Initialisation cancelled."
            << endl;
    }

    // the field has been initialised
    initialised() = true;

    initCoupledBoundaries();
}



// ************************************************************************* //
