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

#include "patchSwitchFvMesh/patchSwitches/trapDoor/trapDoor.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(trapDoor, 0);

    addToRunTimeSelectionTable
    (
        patchSwitch,
        trapDoor,
        dictionary
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool trapDoor::enableCondition()
{
    //Be careful here!
    //The GIBs are disabled. If you compute values on the GIB
    //they will depend on the BC


    //- Check if in first timestep the U field was found in memory
    //  This can happen because the user might want to use a field that is
    //  constructed and stored after mesh (like using a FO)
    if (mesh_.time().restartTimeIndex() == 1)
    {
        if (!mesh_.foundObject<volVectorField>(UName_))
        {
            Info<< "trapDoor "
                 << mesh_.faceZones()[zoneID_].name()
                 << ": velocity field " << UName_
                 << " was not found in memory. "
                 << "Leaving the GIB disabled for this timestep." << endl;
            return false;
        }
    }

    const volVectorField& U =
        mesh_.lookupObject<volVectorField>(UName_);

    const surfaceScalarField& magSf = mesh_.magSf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const scalarField& magSfm = magSf.boundaryField()[masterID_];
    const vectorField& Sfm = Sf.boundaryField()[masterID_];

    const labelList& fcs = mesh_.boundary()[masterID_].faceCells();
    vectorField Umc(U.boundaryField()[masterID_].size(), vector::zero);

    forAll(Umc, cI)
    {
        Umc[cI] = U[fcs[cI]];
    }

    scalar meanUm = 0.0;

    if (mesh_.foundObject<volScalarField>(fractionName_))
    {
        //- Assume only the velocities of the fuel
        const volScalarField& alpha =
            mesh_.lookupObject<volScalarField>(fractionName_);
        const scalarField& alpham = alpha.boundaryField()[masterID_];

        meanUm = gSum
        (
            (Umc*alpham)&Sfm
        )/gSum(magSfm);
    }
    else
    {
        //- Assume single phase (alpha does not exist or user did not specify
        // the phase name correctly)
        meanUm = gSum
        (
            Umc&Sfm
        )/gSum(magSfm);
    }

    if (invertDirection_)
    {
        meanUm = -meanUm;
    }

    Info<< this->type()
         << " for GIB faceZone "
         << mesh_.faceZones()[zoneID_].name()
         << ": open, disableCondition - mean normal fluid velocity: "
         << meanUm << endl;

    return (meanUm>-fVel_);
}

bool trapDoor::disableCondition()
{
    //- Check if in first timestep the U field was found in memory
    //  This can happen because the user might want to use a field that is
    //  constructed and stored after mesh (like using a FO)
    if (mesh_.time().restartTimeIndex() == 1)
    {
        if (!mesh_.foundObject<volScalarField>(pName_))
        {
            Info<< "trapDoor "
                 << mesh_.faceZones()[zoneID_].name()
                 << ": pressure field " << pName_
                 << " was not found in memory. "
                 << "Leaving the GIB enabled for this timestep." << endl;
            return false;
        }
    }

    const volScalarField& p =
        mesh_.lookupObject<volScalarField>(pName_);
    const surfaceScalarField& magSf = mesh_.magSf();

    const scalarField& magSfm = magSf.boundaryField()[masterID_];
    const scalarField& magSfs = magSf.boundaryField()[slaveID_];

    const scalarField& pm = p.boundaryField()[masterID_];
    const scalarField& ps = p.boundaryField()[slaveID_];

    scalar meanPm = gSum
    (
        pm*magSfm
    )/gSum(magSfm);

    scalar meanPs = gSum
    (
        ps*magSfs
    )/gSum(magSfs);

    scalar meanDP = meanPs - meanPm;

    if (invertDirection_)
    {
        meanDP = -meanDP;
    }

    Info<< this->type()
         << " for GIB faceZone "
         << mesh_.faceZones()[zoneID_].name()
         << ": closed, enableCondition - pressureDrop: "
         << meanDP << endl;

    return (meanDP>dp_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

trapDoor::trapDoor
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    patchSwitch(mesh, dict),
    invertDirection_(dict.lookupOrDefault<Switch>("invertDirection", false)),
    fVel_(readScalar(dict.lookup("closeDoorFuelVelocity"))),
    dp_(readScalar(dict.lookup("openDoorPressureDrop"))),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    fractionName_(dict.lookupOrDefault<word>("alpha", "alpha"))
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
