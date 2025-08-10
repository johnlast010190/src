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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "patchSwitchFvMesh/patchSwitches/openPatch/openPatch.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(openPatch, 0);

    addToRunTimeSelectionTable
    (
        patchSwitch,
        openPatch,
        dictionary
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool openPatch::disableCondition()
{
    const volScalarField& field =
        mesh_.lookupObject<volScalarField>(fieldName_);
    const surfaceScalarField& magSf = mesh_.magSf();

    const scalarField& magSfm = magSf.boundaryField()[patchID_];
    const scalarField& fieldpm = field.boundaryField()[patchID_];

    scalar meanFieldpm = gSum
    (
        fieldpm*magSfm
    )/gSum(magSfm);

    Info<< this->type()
         << " for GIB faceZone "
         << mesh_.faceZones()[zoneID_].name()
         << ": closed, enableCondition - field "
         << fieldName_ << " greater than : "
         << meanFieldpm << endl;

    return (meanFieldpm>fieldLim_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

openPatch::openPatch
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    patchSwitch(mesh, dict),
    fieldLim_(readScalar(dict.lookup("triggerValue"))),
    fieldName_(dict.lookupOrDefault<word>("fieldName", "T")),
    side_(dict.lookupOrDefault<word>("side", "master")),
    patchID_(-1)
{
    if (side_ == "master")
    {
        patchID_ = masterID_;
    }
    else if (side_ == "slave")
    {
        patchID_ = slaveID_;
    }
    else
    {
        FatalError
            << "side value of " << side_ << "is not valid." << endl
            << "Only master or slave is allowed"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
