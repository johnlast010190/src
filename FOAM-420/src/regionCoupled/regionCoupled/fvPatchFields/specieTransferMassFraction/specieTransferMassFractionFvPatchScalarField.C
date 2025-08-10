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

#include "specieTransferMassFractionFvPatchScalarField.H"
#include "fvPatchFields/specieTransferVelocity/specieTransferVelocityFvPatchVectorField.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "basicThermo/basicThermo.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum
    <
        specieTransferMassFractionFvPatchScalarField::property,
        4
    >::names[] =
    {
        "massFraction",
        "moleFraction",
        "molarConcentration",
        "partialPressure"
    };
}

const Foam::NamedEnum
<
    Foam::specieTransferMassFractionFvPatchScalarField::property,
    4
> Foam::specieTransferMassFractionFvPatchScalarField::propertyNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specieTransferMassFractionFvPatchScalarField::
specieTransferMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    UName_("U"),
    phiYp_(p.size(), 0),
    timeIndex_(-1),
    clipMassFractions_(true),
    c_(0),
    property_(massFraction)
{}


Foam::specieTransferMassFractionFvPatchScalarField::
specieTransferMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF, dict, false),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiYp_(p.size(), 0),
    timeIndex_(-1),
    clipMassFractions_(dict.lookupOrDefault<Switch>("clipMassFraction", true)),
    c_(dict.lookupOrDefault<scalar>("c", scalar(0))),
    property_
    (
        c_ == scalar(0)
      ? massFraction
      : propertyNames_.read(dict.lookup("property"))
    )
{
    forceAssign(scalarField("value", dict, p.size()));

    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::specieTransferMassFractionFvPatchScalarField::
specieTransferMassFractionFvPatchScalarField
(
    const specieTransferMassFractionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_),
    phiYp_(p.size(), 0),
    timeIndex_(-1),
    clipMassFractions_(ptf.clipMassFractions_),
    c_(ptf.c_),
    property_(ptf.property_)
{}


Foam::specieTransferMassFractionFvPatchScalarField::
specieTransferMassFractionFvPatchScalarField
(
    const specieTransferMassFractionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_),
    phiYp_(ptf.size(), 0),
    timeIndex_(-1),
    clipMassFractions_(ptf.clipMassFractions_),
    c_(ptf.c_),
    property_(ptf.property_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::scalarField&
Foam::specieTransferMassFractionFvPatchScalarField::phiYp() const
{
    if (timeIndex_ != this->db().time().timeIndex())
    {
        timeIndex_ = this->db().time().timeIndex();

        phiYp_ = calcPhiYp();
    }

    return phiYp_;
}


void Foam::specieTransferMassFractionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    mixedFvPatchField<scalar>::evaluate(commsType);

    if (clipMassFractions_)
    {
        // Additional clipping of the value fractions to [0,1]
        forAll(patch(), facei)
        {
            if ((*this)[facei] > 1)
            {
                (*this)[facei] = 1;
            }
            else if ((*this)[facei] < 0)
            {
                (*this)[facei] = 0;
            }
        }
    }
}


void Foam::specieTransferMassFractionFvPatchScalarField::updateCoeffs()
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

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );
    const scalarField muEffp(turbModel.muEff(patch().index()));
    const scalarField AMuEffp(patch().magSf()*muEffp);

    // Set the gradient and value so that the transport and diffusion combined
    // result in the desired specie flux
    valueFraction() = phip/(phip - patch().deltaCoeffs()*AMuEffp);
    refValue() = *this;
    refGrad() = phip*(*this - phiYp()/uPhip)/AMuEffp;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::specieTransferMassFractionFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    writeEntryIfDifferent<scalar>(os, "c", scalar(0), c_);
    os.writeEntry("property", propertyNames_[property_]);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<Switch>
    (
        os,
        "clipMassFractions",
        true,
        clipMassFractions_
    );
}


// ************************************************************************* //
