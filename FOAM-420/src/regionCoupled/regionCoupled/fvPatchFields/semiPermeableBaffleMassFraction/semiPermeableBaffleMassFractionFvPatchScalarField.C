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
    (c) 2017-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvPatchFields/semiPermeableBaffleMassFraction/semiPermeableBaffleMassFractionFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "basicThermo/basicThermo.H"
#include "mappedPatches/mappedPolyPatch/mappedPatchBase.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/mapped/mappedFvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::semiPermeableBaffleMassFractionFvPatchScalarField::
semiPermeableBaffleMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    specieTransferMassFractionFvPatchScalarField(p, iF),
    matTables_(nullptr),
    nbrMatTables_(nullptr),
    turb_(nullptr),
    nbrTurb_(nullptr),
    diffusionLimited_(true)
{}


Foam::semiPermeableBaffleMassFractionFvPatchScalarField::
semiPermeableBaffleMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    specieTransferMassFractionFvPatchScalarField(p, iF, dict),
    matTables_(nullptr),
    nbrMatTables_(nullptr),
    turb_(nullptr),
    nbrTurb_(nullptr),
    diffusionLimited_(dict.lookupOrDefault<Switch>("diffusionLimited", true))
{}


Foam::semiPermeableBaffleMassFractionFvPatchScalarField::
semiPermeableBaffleMassFractionFvPatchScalarField
(
    const semiPermeableBaffleMassFractionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    specieTransferMassFractionFvPatchScalarField(ptf, p, iF, mapper),
    matTables_(nullptr),
    nbrMatTables_(nullptr),
    turb_(nullptr),
    nbrTurb_(nullptr),
    diffusionLimited_(ptf.diffusionLimited_)
{}


Foam::semiPermeableBaffleMassFractionFvPatchScalarField::
semiPermeableBaffleMassFractionFvPatchScalarField
(
    const semiPermeableBaffleMassFractionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    specieTransferMassFractionFvPatchScalarField(ptf, iF),
    matTables_(nullptr),
    nbrMatTables_(nullptr),
    turb_(nullptr),
    nbrTurb_(nullptr),
    diffusionLimited_(ptf.diffusionLimited_)
{}


// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//

void
Foam::semiPermeableBaffleMassFractionFvPatchScalarField::lookupProperties() const
{
    const mappedFvPatch& pMap = refCast<const mappedFvPatch>(patch());
    const polyMesh& nbrMesh = pMap.nbrMesh();

    if (!turb_)
    {
        turb_ =
            &db().lookupObject<cmpTurbModel>
            (
                IOobject::groupName
                (
                    cmpTurbModel::propertiesName,
                    internalField().group()
                )
            );
    }

    if (!nbrTurb_)
    {
        nbrTurb_ =
            &nbrMesh.thisDb().lookupObject<cmpTurbModel>
            (
                IOobject::groupName
                (
                    cmpTurbModel::propertiesName,
                    internalField().group()
                )
            );
    }

    if (!nbrMatTables_)
    {
        nbrMatTables_ =
            &nbrMesh.subRegistry("materialModels").lookupObject<materialTables>
            (
                "materialTables"
            );
    }

    if (!matTables_)
    {
        matTables_ =
            &db().subRegistry("materialModels").lookupObject<materialTables>
            (
                "materialTables"
            );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::semiPermeableBaffleMassFractionFvPatchScalarField::calcPhiYp() const
{
    const label nFaces = patch().size();
    if (c_ == scalar(0))
    {
        return tmp<scalarField>(new scalarField(nFaces, Zero));
    }

    const label patchi = patch().index();
    const word& YName = internalField().name();

    // Coupling information
    const mappedFvPatch& pMap = refCast<const mappedFvPatch>(patch());
    const label nbrPatchi = pMap.nbrPatchID();
    const polyMesh& nbrMesh = pMap.nbrMesh();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[nbrPatchi];
    const word phaseName(internalField().group());

    const scalarField Yc(patchInternalField());
    const scalarField nbrYc
    (
        pMap.interpolate
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>(YName)
           .patchInternalField(),
            scalarField(nFaces, Zero)
        )
    );

    // Get the specie molecular weight, if needed
    scalar Wi = NaN;
    if (property_ != massFraction)
    {
        Wi = tables()(WModel::typeName, phaseName, YName)[0];
    }

    // Get the mixture molecular weights, if needed
    tmp<scalarField> tW, tNbrW;
    if (property_ == moleFraction || property_ == partialPressure)
    {
        tW = tables()(WModel::typeName).boundaryField()[patchi];
        tNbrW =
            pMap.interpolate
            (
                nbrTables()(WModel::typeName).boundaryField()[nbrPatchi],
                scalarField(nFaces, Zero)
            );
    }

    // Construct coefficients that convert mass fraction to the property that
    // drives the transfer
    scalarField k(nFaces, 1), nbrK(nFaces, 1);
    switch(property_)
    {
        case massFraction:
        {
            break;
        }

        case moleFraction:
        {
            k *= tW/Wi;
            nbrK *= tNbrW/Wi;
            break;
        }

        case molarConcentration:
        {
            k *= tables()(rhoModel::typeName).boundaryField()[patchi]/Wi;
            tmp<scalarField> nbrRhop =
                pMap.interpolate
                (
                    nbrTables()(rhoModel::typeName).boundaryField()[nbrPatchi],
                    scalarField(nFaces, Zero)
                );
            nbrK *= nbrRhop/Wi;
            break;
        }

        case partialPressure:
        {
            const scalarField p
            (
                db().lookupObject<refScalarField>("pRef").patchField(patchi)
            );
            const scalarField pNei
            (
                nbrMesh.thisDb().lookupObject<refScalarField>("pRef")
               .patchField(nbrPatchi)
            );

            k *= p*tW/Wi;
            tmp<scalarField> nbrPp =
                pMap.interpolate(pNei, scalarField(nFaces, Zero));
            nbrK *= nbrPp*tNbrW/Wi;
            break;
        }
    }

    if (diffusionLimited_)
    {
        const scalarField alphaEffDeltap
        (
            turb().kappaEff(patchi)*patch().deltaCoeffs()
           /tables()(CpModel::typeName).boundaryField()[patchi]
        );

        const scalarField nbrAlphaEffDeltap
        (
            pMap.interpolate
            (
                nbrTurb().kappaEff(nbrPatchi)*nbrPatch.deltaCoeffs()
               /nbrTables()(CpModel::typeName).boundaryField()[nbrPatchi],
               scalarField(nFaces, Zero)
            )
        );
        return
            (patch().magSf()*(k*Yc - nbrK*nbrYc))
           /(1/c_ + k/alphaEffDeltap + nbrK/nbrAlphaEffDeltap);
    }

    return patch().magSf()*c_*(k*Yc - nbrK*nbrYc);
}


void Foam::semiPermeableBaffleMassFractionFvPatchScalarField::write
(
    Ostream& os
) const
{
    specieTransferMassFractionFvPatchScalarField::write(os);
    writeEntryIfDifferent<Switch>
    (
        os,
        "diffusionLimited",
        true,
        diffusionLimited_
    );
}



// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        semiPermeableBaffleMassFractionFvPatchScalarField
    );
}

// ************************************************************************* //
