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
    (c) 2011 OpenFOAM Foundation
    (c) 2010-2018 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/radiant/radiantFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "fvMesh/fvPatches/derived/mapped/mappedFvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiantFvPatchScalarField::
radiantFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    transmissivity_(p.size(), 0.0),
    radiationFlux_(0.0)
{}

Foam::radiantFvPatchScalarField::
radiantFvPatchScalarField
(
    const radiantFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    transmissivity_(mapper(ptf.transmissivity_)),
    radiationFlux_(ptf.radiationFlux_)
{}

Foam::radiantFvPatchScalarField::
radiantFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF, dict),
    transmissivity_(p.size(), 1.0),
    radiationFlux_(dict.lookupOrDefault<scalar>("radiationFlux", 0.0))
{
    if (dict.found("transmissivity"))
    {
        transmissivity_ = scalarField("transmissivity", dict, p.size());
    }
}

Foam::radiantFvPatchScalarField::
radiantFvPatchScalarField
(
    const radiantFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf),
    transmissivity_(ptf.transmissivity_),
    radiationFlux_(ptf.radiationFlux_)
{}


Foam::radiantFvPatchScalarField::
radiantFvPatchScalarField
(
    const radiantFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    transmissivity_(ptf.transmissivity_),
    radiationFlux_(ptf.radiationFlux_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::radiantFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    refValue() = radiationFlux_ * transmissivity_;

    if (isA<mappedFvPatch>(this->patch()))
    {
        // Since we're inside initEvaluate/evaluate there might be processor
        // comms underway. Change the tag we use.
        int oldTag = UPstream::msgType();
        UPstream::msgType() = oldTag+1;

        const mappedFvPatch& mpp = refCast<const mappedFvPatch>
        (
            patch()
        );

        const fvPatch& nbrPatch = mpp.nbrPatch();

        //we need to check if this field exists on the other side
        //of the mapped patch otherwise the solver will fail for
        //multi-region cases
        const polyMesh& nbrMesh = mpp.nbrMesh();
        const polyMesh& mesh = patch().boundaryMesh().mesh();

        if
        (
            nbrMesh == mesh
            ||
            (
                nbrMesh != mesh &&
                nbrMesh.foundObject<volScalarField>(internalField().name())
            )
        )
        {
            scalarField nbrFlux =
                nbrPatch.lookupPatchField<volScalarField, scalar>
                (this->internalField().name());
            const scalarField defaultValues(patch().size(), 0.0);
            nbrFlux = mpp.interpolate
            (
                nbrFlux,
                defaultValues
            );

            refValue() = transmissivity_*nbrFlux;
        }
        else
        {
            // no sunshine in neighbour region
            refValue() = 0.0;
        }
    }


    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::radiantFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    if (phiName_ != "phi")
    {
        os.writeEntry("phi", phiName_);
    }
    this->refValue().writeEntry("inletValue", os);
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::radiantFvPatchScalarField::operator=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchField<scalar>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        radiantFvPatchScalarField
    );
}

// ************************************************************************* //
