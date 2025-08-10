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
    (c) 2019-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvPatchFields/specieTransferVelocity/specieTransferVelocityFvPatchVectorField.H"
#include "fvPatchFields/specieTransferMassFraction/specieTransferMassFractionFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "basicThermo/basicThermo.H"
#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "rhoReactionThermo/rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specieTransferVelocityFvPatchVectorField::
specieTransferVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueInletOutletFvPatchVectorField(p, iF),
    rhoName_("rho")
{}


Foam::specieTransferVelocityFvPatchVectorField::
specieTransferVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueInletOutletFvPatchVectorField(p, iF, dict, false),
    rhoName_(dict.lookupOrDefault<word>("rhoName", "rho"))
{
    if (dict.found("value"))
    {
        forceAssign(vectorField("value", dict, p.size()));
    }
    else
    {
        specieTransferVelocityFvPatchVectorField::updateCoeffs();
    }
}


Foam::specieTransferVelocityFvPatchVectorField::
specieTransferVelocityFvPatchVectorField
(
    const specieTransferVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueInletOutletFvPatchVectorField(ptf, p, iF, mapper),
    rhoName_(ptf.rhoName_)
{}


Foam::specieTransferVelocityFvPatchVectorField::
specieTransferVelocityFvPatchVectorField
(
    const specieTransferVelocityFvPatchVectorField& rwvpvf
)
:
    fixedValueInletOutletFvPatchVectorField(rwvpvf),
    rhoName_(rwvpvf.rhoName_)
{}


Foam::specieTransferVelocityFvPatchVectorField::
specieTransferVelocityFvPatchVectorField
(
    const specieTransferVelocityFvPatchVectorField& rwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueInletOutletFvPatchVectorField(rwvpvf, iF),
    rhoName_(rwvpvf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
const Foam::tmp<Foam::scalarField>
Foam::specieTransferVelocityFvPatchVectorField::phip() const
{
    typedef specieTransferMassFractionFvPatchScalarField YBCType;

    const basicSpecieMixture& composition =
        db().lookupObject<rhoReactionThermo>
        (
            basicThermo::matDictName
        ).composition();
    const PtrList<volScalarField>& Y = composition.Y();

    // Sum up the phiYp-s from all the species
    tmp<scalarField> tPhip(new scalarField(this->size(), 0));
    scalarField& phip = tPhip.ref();
    forAll(Y, i)
    {
        const fvPatchScalarField& Yp = Y[i].boundaryField()[patch().index()];

        if (!isA<YBCType>(Yp))
        {
            FatalErrorInFunction
                << "The mass-fraction condition on patch " << patch().name()
                << " is not of type " << YBCType::typeName << "."
                << exit(FatalError);
        }

        phip += refCast<const YBCType>(Yp).phiYp();
    }

    return tPhip;
}

void Foam::specieTransferVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueInletOutletFvPatchVectorField::autoMap(m);
}


void Foam::specieTransferVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueInletOutletFvPatchVectorField::rmap(ptf, addr);
}


void Foam::specieTransferVelocityFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueInletOutletFvPatchVectorField::autoMapGIB(mapper);
}


void Foam::specieTransferVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the density
    const scalarField& rhop =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    // Set the normal component of the velocity to match the computed flux
    const vectorField nf(patch().nf());
    const tensorField Tau(tensor::I - sqr(nf));
    this->forceAssign((Tau & *this) + nf*phip()/(rhop*patch().magSf()));

    fixedValueInletOutletFvPatchVectorField::updateCoeffs();
}


void Foam::specieTransferVelocityFvPatchVectorField::write(Ostream& os) const
{
    fixedValueInletOutletFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        specieTransferVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
