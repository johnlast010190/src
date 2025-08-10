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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "wallContactTime/wallContactTime.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallContactTime, 0);
    addToRunTimeSelectionTable(functionObject, wallContactTime, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallContactTime::wallContactTime
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    patchNames_(0),
    patchIDs_(0),
    fieldNames_(0)
{
    read(dict);
    createFields(phaseContactTimes_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallContactTime::~wallContactTime()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::functionObjects::wallContactTime::createFields
(
    PtrList<volScalarField>& flds
) const
{
    Log << type() << " " << name() <<  " createFields:" << nl;
    flds.resize(fieldNames_.size());
    forAll(fieldNames_, i)
    {
        const word& fldName = word(fieldNames_[i] + "ContactTime");
        const volScalarField& alpha =
            mesh_.lookupObject<volScalarField>(fieldNames_[i]);

        flds.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fldName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(fldName,alpha.dimensions()*dimTime,0),
                fixedValueFvPatchScalarField::typeName
            )
        );
    }
}

bool Foam::functionObjects::wallContactTime::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    //calculate
    calculate();

    Log <<  "Calculated wall contact time " << endl;
    Log << endl;

    return true;
}


void Foam::functionObjects::wallContactTime::calculate()
{
    // update wall contact time
    forAll(fieldNames_, i)
    {
        const volScalarField& alpha =
            mesh_.lookupObject<volScalarField>(fieldNames_[i]);

        scalar deltaT = mesh_.time().deltaT().value();

        forAll(patchIDs_, pI)
        {
            if (patchIDs_[pI] != -1)
            {
                fvPatchField<scalar>& pfld =
                    phaseContactTimes_[i].boundaryFieldRef()[patchIDs_[pI]];
                pfld.forceAssign(pfld + alpha.boundaryField()[patchIDs_[pI]]*deltaT);
            }
        }
    }
}

bool Foam::functionObjects::wallContactTime::write()
{
    Log << type() << " " << name() <<  " write:" << nl;

    //calculate statistics
    forAll(fieldNames_, i)
    {
        forAll(patchIDs_, pI)
        {
            if (patchIDs_[pI] != -1)
            {
                const fvPatchField<scalar>& pfld =
                    phaseContactTimes_[i].boundaryField()[patchIDs_[pI]];

                scalar fldMin(min(pfld));
                scalar fldMax(max(pfld));
                scalar fldAve(0);
                scalar totalArea(0);

                const scalarField& area
                    = mesh_.magSf().boundaryField()[patchIDs_[pI]];

                forAll(pfld, fI)
                {
                    fldAve += pfld[fI] * area[fI];
                    totalArea += area[fI];
                }

                reduce(fldMin, minOp<scalar>());
                reduce(fldMax, maxOp<scalar>());
                reduce(fldAve, sumOp<scalar>());
                reduce(totalArea, sumOp<scalar>());

                fldAve /= totalArea;

                Log << "min/max/mean contact time for phase" << fieldNames_[i]
                    << " on patch " << pfld.patch().name()
                    << " : " << fldMin << " , " << fldMax << " , " << fldAve << endl;
            }
        }
        Log << endl;
    }

    return true;
}


bool Foam::functionObjects::wallContactTime::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name() <<  " read:" << nl;

    patchNames_ = wordList(dict.lookup("patches"));
    patchIDs_.setSize(patchNames_.size());

    forAll(patchNames_, pI)
    {
        patchIDs_[pI] = mesh_.boundaryMesh().findPatchID(patchNames_[pI]);

        if (patchIDs_[pI] == -1)
        {
            WarningInFunction
                << "Could not find patch with name "
                << patchNames_[pI] << " in boundary list." << nl
                << "Entry will be ignored by wallContactTime:"<< name() << "."
                << endl;
        }
    }

    fieldNames_ = wordList(dict.lookup("fields"));

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
