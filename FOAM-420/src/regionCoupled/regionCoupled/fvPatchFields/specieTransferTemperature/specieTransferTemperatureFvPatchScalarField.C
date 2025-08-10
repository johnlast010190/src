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

#include "fvPatchFields/specieTransferTemperature/specieTransferTemperatureFvPatchScalarField.H"
#include "fvPatchFields/specieTransferMassFraction/specieTransferMassFractionFvPatchScalarField.H"
#include "fvPatchFields/specieTransferVelocity/specieTransferVelocityFvPatchVectorField.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "rhoReactionThermo/rhoReactionThermo.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::specieTransferTemperatureFvPatchScalarField::setThermos() const
{
    if (!thermoPtr_)
    {
        thermoPtr_ =
        (
            &db().lookupObject<basicThermo>
            (
                IOobject::groupName
                (
                    basicThermo::dictName,
                    internalField().group()
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specieTransferTemperatureFvPatchScalarField::
specieTransferTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedEnergyCalculatedTemperatureFvPatchScalarField(p, iF),
    phiName_("phi"),
    UName_("U"),
    refTempGrad_(p.size(), 0),
    refTemp_(p.size(), 0),
    thermoPtr_(nullptr)
{}


Foam::specieTransferTemperatureFvPatchScalarField::
specieTransferTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool readValue
)
:
    mixedEnergyCalculatedTemperatureFvPatchScalarField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    refTempGrad_(p.size(), 0),
    refTemp_(p.size(), 0),
    thermoPtr_(nullptr)
{
    if (readValue)
    {
        forceAssign(scalarField("value", dict, p.size()));
    }
}


Foam::specieTransferTemperatureFvPatchScalarField::
specieTransferTemperatureFvPatchScalarField
(
    const specieTransferTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedEnergyCalculatedTemperatureFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_),
    refTempGrad_(p.size(), 0),
    refTemp_(p.size(), 0),
    thermoPtr_(nullptr)
{}


Foam::specieTransferTemperatureFvPatchScalarField::
specieTransferTemperatureFvPatchScalarField
(
    const specieTransferTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedEnergyCalculatedTemperatureFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_),
    refTempGrad_(ptf.size(), 0),
    refTemp_(ptf.size(), 0),
    thermoPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::tmp<Foam::scalarField>
Foam::specieTransferTemperatureFvPatchScalarField::phiHep() const
{
    typedef specieTransferMassFractionFvPatchScalarField YBCType;

    const word phaseName(internalField().group());
    const materialTables& matTables =
        db().subRegistry("materialModels").lookupObject<materialTables>
        (
            "materialTables"
        );
    const basicSpecieMixture& composition =
        db().lookupObject<rhoReactionThermo>
        (
            basicThermo::matDictName
        ).composition();
    const PtrList<volScalarField>& Y = composition.Y();

    // Sum up the phiHep from all the species
    tmp<scalarField> tPhiHep(new scalarField(size(), 0));
    scalarField& phiHep = tPhiHep.ref();
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

        // Note: Temperature and pressure will be already calculated based
        // on field boundary value of temperature and pressure.
        phiHep +=
            refCast<const YBCType>(Yp).phiYp()
           *matTables
            (
                HEModel::typeName,
                phaseName,
                Y[i].name()
            ).boundaryField()[patch().index()];
    }

    return tPhiHep;
}


void Foam::specieTransferTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the fluxes
    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);
    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);
    tmp<scalarField> uPhip =
        refCast<const specieTransferVelocityFvPatchVectorField>(Up).phip();

    const compressible::turbulenceModel& turbModel =
        db().lookupObject<compressible::turbulenceModel>
        (
            IOobject::groupName
            (
                compressible::turbulenceModel::propertiesName,
                internalField().group()
            )
        );
    setThermos();

    const label patchi = patch().index();
    // Get the diffusivity
    const scalarField Cp(thermoPtr_->Cp().ref().boundaryField()[patchi]);
    const scalarField AAlphaEffp(patch().magSf()*turbModel.kappaEff(patchi));

    // Get the current energy to linearise around
    const scalarField& hep = thermoPtr_->he().boundaryField()[patchi];

    if (ishe())
    {
        heValueFraction() = (phip*Cp)/(phip - patch().deltaCoeffs()*AAlphaEffp);
        heRefValue() = hep;
        refTemp_ = hep/Cp;
        refTempGrad_ = phip*(hep - phiHep()/uPhip)/AAlphaEffp;
        heRefGrad() = refTempGrad_*Cp;

        mixedEnergyCalculatedTemperatureFvPatchScalarField::updateCoeffs();
    }
    else
    {
        fvPatchField<scalar>::updateCoeffs();
    }
}


Foam::tmp<Foam::scalarField>
Foam::specieTransferTemperatureFvPatchScalarField::snGrad() const
{
    if (ishe())
    {
        mixedEnergyCalculatedTemperatureFvPatchScalarField::snGrad();
    }
    const scalarField Cp
    (
        thermoPtr_->Cp().ref().boundaryField()[patch().index()].internalField()
    );

    return
        heValueFraction()
       *(refTemp_ - (patchInternalField()/Cp))
       *patch().deltaCoeffs()
      + (1.0 - heValueFraction())*refTempGrad_;
}


void Foam::specieTransferTemperatureFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (ishe())
    {
        mixedEnergyCalculatedTemperatureFvPatchScalarField::evaluate
        (
            commsType
        );
    }
}


Foam::tmp<Foam::scalarField>
Foam::specieTransferTemperatureFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    if (ishe())
    {
        mixedEnergyCalculatedTemperatureFvPatchScalarField::
        valueInternalCoeffs(w);
    }
    return (1.0 - heValueFraction());
}


Foam::tmp<Foam::scalarField>
Foam::specieTransferTemperatureFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    if (ishe())
    {
        mixedEnergyCalculatedTemperatureFvPatchScalarField::
        valueBoundaryCoeffs(w);
    }

    return
        heValueFraction()*refTemp_
      + (1.0 - heValueFraction())*refTempGrad_/patch().deltaCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::specieTransferTemperatureFvPatchScalarField::
gradientInternalCoeffs() const
{
    if (ishe())
    {
        mixedEnergyCalculatedTemperatureFvPatchScalarField::
        gradientInternalCoeffs();
    }
    return -heValueFraction()*patch().deltaCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::specieTransferTemperatureFvPatchScalarField::
gradientBoundaryCoeffs() const
{
    if (ishe())
    {
        mixedEnergyCalculatedTemperatureFvPatchScalarField::
        gradientBoundaryCoeffs();
    }
    return
        heValueFraction()*patch().deltaCoeffs()*refTemp_
      + (1.0 - heValueFraction())*refTempGrad_;
}


void Foam::specieTransferTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        specieTransferTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
