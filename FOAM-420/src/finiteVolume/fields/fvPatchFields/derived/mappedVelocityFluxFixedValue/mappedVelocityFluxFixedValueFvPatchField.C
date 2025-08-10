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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/mappedVelocityFluxFixedValue/mappedVelocityFluxFixedValueFvPatchField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "mappedPatches/mappedPolyPatch/mappedPatchBase.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedVelocityFluxFixedValueFvPatchField::
mappedVelocityFluxFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi")
{}


Foam::mappedVelocityFluxFixedValueFvPatchField::
mappedVelocityFluxFixedValueFvPatchField
(
    const mappedVelocityFluxFixedValueFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "Patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << " for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }
}


Foam::mappedVelocityFluxFixedValueFvPatchField::
mappedVelocityFluxFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "Patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << " for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        this->patch().patch()
    );
    if (mpp.mode() == mappedPolyPatch::NEARESICELL)
    {
        FatalErrorInFunction
            << "Patch " << p.name()
            << " of type '" << p.type()
            << "' can not be used in 'nearestCell' mode"
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }
}


Foam::mappedVelocityFluxFixedValueFvPatchField::
mappedVelocityFluxFixedValueFvPatchField
(
    const mappedVelocityFluxFixedValueFvPatchField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    phiName_(ptf.phiName_)
{}


Foam::mappedVelocityFluxFixedValueFvPatchField::
mappedVelocityFluxFixedValueFvPatchField
(
    const mappedVelocityFluxFixedValueFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedVelocityFluxFixedValueFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        mappedVelocityFluxFixedValueFvPatchField::patch().patch()
    );
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    const word& fieldName = internalField().name();
    const volVectorField& UField =
        nbrMesh.lookupObject<volVectorField>(fieldName);

    const surfaceScalarField& phiField =
        nbrMesh.lookupObject<surfaceScalarField>(phiName_);

    vectorField newUValues;
    scalarField newPhiValues;

    switch (mpp.mode())
    {
        case mappedPolyPatch::NEARESIFACE:
        {
            vectorField allUValues(nbrMesh.nFaces(), Zero);
            scalarField allPhiValues(nbrMesh.nFaces(), 0.0);

            forAll(UField.boundaryField(), patchi)
            {
                const fvPatchVectorField& Upf = UField.boundaryField()[patchi];
                const scalarField& phipf = phiField.boundaryField()[patchi];

                label faceStart = Upf.patch().start();

                forAll(Upf, facei)
                {
                    allUValues[faceStart + facei] = Upf[facei];
                    allPhiValues[faceStart + facei] = phipf[facei];
                }
            }

            mpp.distribute(allUValues);
            newUValues.transfer(allUValues);

            mpp.distribute(allPhiValues);
            newPhiValues.transfer(allPhiValues);

            break;
        }
        case mappedPolyPatch::NEARESIPATCHFACE:
        case mappedPolyPatch::NEARESIPATCHFACEAMI:
        {
            const label nbrPatchID =
                nbrMesh.boundaryMesh().findPatchID(mpp.samplePatch());

            newUValues = UField.boundaryField()[nbrPatchID];
            mpp.distribute(newUValues);

            newPhiValues = phiField.boundaryField()[nbrPatchID];
            mpp.distribute(newPhiValues);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "patch can only be used in NEARESIPATCHFACE, "
                << "NEARESIPATCHFACEAMI or NEARESIFACE mode" << nl
                << abort(FatalError);
        }
    }

    forceAssign(newUValues);
    const_cast<surfaceScalarField&>
    (
        phiField
    ).boundaryFieldRef()[patch().index()].forceAssign(newPhiValues);

    // Restore tag
    UPstream::msgType() = oldTag;

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::mappedVelocityFluxFixedValueFvPatchField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        mappedVelocityFluxFixedValueFvPatchField
    );
}


// ************************************************************************* //
