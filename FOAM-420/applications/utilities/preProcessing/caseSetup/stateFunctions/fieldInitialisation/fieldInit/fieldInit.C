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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldInit.H"
#include "fields/volFields/volFields.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldInit, 0);
    defineRunTimeSelectionTable(fieldInit, initMethod);
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldInit::fieldInit
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDefinition,
    const word& fN
)
:
mesh_(mesh),
localDb_(localDb),
fieldName_(fN),
fieldDict_(fieldDefinition),
initialized_(false)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::dictionary& Foam::fieldInit::initDict() const
{
    return fieldDict_.subDict("initialisation");
}

void Foam::fieldInit::initMsg() const
{
    Info<< "   - " <<  name() << " : " << type() << endl;
}

void Foam::fieldInit::initTurbBoundaries() const
{
    volScalarField& f
        = const_cast<volScalarField&>
        (localDb().lookupObject<volScalarField>(fieldName_));

    forAll(mesh_.boundary(), pI)
    {
        if
        (
            isA<calculatedFvPatchField<scalar>>
            (
                f.boundaryField()[pI]
            )
            ||
            isA<mixedFvPatchField<scalar>>
            (
                f.boundaryField()[pI]
            )
        )
        {
            f.boundaryFieldRef()[pI].forceAssign
            (
                f.boundaryField()[pI].patchInternalField()
            );
        }
    }
}

void Foam::fieldInit::initCoupledBoundaries() const
{

    //apply zero gradient to coupled boundaries
    //update coupled patches - prevents turb properties from being zero
    //after parallel initialisation

    if (localDb().foundObject<volSphericalTensorField>(fieldName_))
    {
        volSphericalTensorField& f
            = const_cast<volSphericalTensorField&>
            (localDb().lookupObject<volSphericalTensorField>(fieldName_));

        forAll(mesh_.boundary(), pI)
        {
            if (mesh_.boundary()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }
            else if (f.boundaryField()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }

        }

    }
    else if (localDb().foundObject<volSymmTensorField>(fieldName_))
    {
        volSymmTensorField& f
            = const_cast<volSymmTensorField&>
            (localDb().lookupObject<volSymmTensorField>(fieldName_));

        forAll(mesh_.boundary(), pI)
        {
            if (mesh_.boundary()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }
            else if (f.boundaryField()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }
        }

    }
    else if (localDb().foundObject<volTensorField>(fieldName_))
    {
        volTensorField& f
            = const_cast<volTensorField&>
            (localDb().lookupObject<volTensorField>(fieldName_));

        forAll(mesh_.boundary(), pI)
        {
            if (mesh_.boundary()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }
            else if (f.boundaryField()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }
        }

    }
    else if (localDb().foundObject<volVectorField>(fieldName_))
    {
        volVectorField& f
            = const_cast<volVectorField&>
            (localDb().lookupObject<volVectorField>(fieldName_));

        forAll(mesh_.boundary(), pI)
        {
            if (mesh_.boundary()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }
            else if (f.boundaryField()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }
        }

    }
    else if (localDb().foundObject<volScalarField>(fieldName_))
    {
        volScalarField& f
            = const_cast<volScalarField&>
            (localDb().lookupObject<volScalarField>(fieldName_));

        forAll(mesh_.boundary(), pI)
        {
            if (mesh_.boundary()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }
            else if (f.boundaryField()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }
        }
    }
    else if (mesh().foundObject<pointScalarField>(fieldName_))
    {
        pointScalarField& f
            = const_cast<pointScalarField&>
            (mesh().lookupObject<pointScalarField>(fieldName_));

        forAll(mesh_.boundary(), pI)
        {
            if (mesh_.boundary()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }
            else if (f.boundaryField()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }
        }
    }
    else if (mesh().foundObject<pointVectorField>(fieldName_))
    {
        pointVectorField& f
            = const_cast<pointVectorField&>
            (mesh().lookupObject<pointVectorField>(fieldName_));

        forAll(mesh_.boundary(), pI)
        {
            if (mesh_.boundary()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );

            }
            else if (f.boundaryField()[pI].coupled())
            {
                f.boundaryFieldRef()[pI].forceAssign(
                    f.boundaryField()[pI].patchInternalField()
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Trying to initialise unsupported field type."
            << exit(FatalError);
    }

}


// ************************************************************************* //
