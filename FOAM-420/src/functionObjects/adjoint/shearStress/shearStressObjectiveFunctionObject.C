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

#include "shearStress/shearStressObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(shearStressObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        shearStressObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::shearStressObjectiveFunctionObject::
shearStressObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    ssSource_(objectiveDict.lookupOrDefault<Switch>("meanShearStress", false)),
    ssMeanPtr_(nullptr)
{
    createFiles(useAdjointFileFormat);

    if (ssSource_)
    {
        ssMeanPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "shearStressMean",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }
}


Foam::functionObjects::shearStressObjectiveFunctionObject::
shearStressObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    shearStressObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::shearStressObjectiveFunctionObject::
~shearStressObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::shearStressObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    return true;
}


bool
Foam::functionObjects::shearStressObjectiveFunctionObject::
execute()
{
    tmp<volScalarField> muEff = objectiveFunctionObject::muEff();
    const volScalarField::Boundary& muEffw
        = muEff().boundaryField();
    tmp<volVectorField> V = U();

    scalar totalObjective = 0;

    forAll(objectivePatch_, pI)
    {
        if (objectivePatch_[pI] && isA<wallFvPatch>(mesh_.boundary()[pI]))
        {
            const scalarField muEffwp = muEffw[pI];
            const fvPatchVectorField& Uw = V().boundaryField()[pI];
            const vectorField& Uw2 = V().boundaryField()[pI];
            vectorField Uint( Uw.patchInternalField() );

            const vectorField nf( mesh_.boundary()[pI].nf() );
            scalarField Uw2n( Uw2 & nf );
            scalarField Uintn( Uint & nf );
            vectorField Uw2t( Uw2 - (Uw2n * nf) );
            vectorField Uintt( Uint - (Uintn * nf) );
            const scalarField& deltaC = mesh_.boundary()[pI].deltaCoeffs();
            vectorField dUdn( muEffwp*((Uw2t - Uintt)*deltaC) );
            scalarField dUdn_2( 0.5*mag(dUdn)*mag(dUdn) );
            totalObjective += sum(dUdn_2);
        }
    }

    reduce(totalObjective, sumOp<scalar>());
    objectiveValue_ = totalObjective;

    Info<< type() << " " << name() << " execute:" << nl
        << "Shear stress = " << objectiveValue_ << " [N]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::shearStressObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
