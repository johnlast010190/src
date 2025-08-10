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
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "relativeHumidityInletOutletFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "basicThermo/basicThermo.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeHumidityInletOutletFvPatchScalarField::
relativeHumidityInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    relativeHumidity_(p.size(), 0.0),
    phiName_("phi"),
    TName_("T"),
    PName_("p"),
    rhoName_("rho"),
    Mvap_(18.02),
    Mmix_(28.96)
{}

Foam::relativeHumidityInletOutletFvPatchScalarField::
relativeHumidityInletOutletFvPatchScalarField
(
    const relativeHumidityInletOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    relativeHumidity_(mapper(ptf.relativeHumidity_)),
    phiName_(ptf.phiName_),
    TName_(ptf.TName_),
    PName_(ptf.PName_),
    rhoName_(ptf.rhoName_),
    Mvap_(ptf.Mvap_),
    Mmix_(ptf.Mmix_)
{}

Foam::relativeHumidityInletOutletFvPatchScalarField::
relativeHumidityInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF, dict, false),
    relativeHumidity_("relativeHumidity", dict, p.size()),
    phiName_(dict.lookupOrDefault<word>("phi","phi")),
    TName_(dict.lookupOrDefault<word>("T","T")),
    PName_(dict.lookupOrDefault<word>("p","p")),
    rhoName_(dict.lookupOrDefault<word>("rho","rho")),
    Mvap_(dict.lookupOrDefault<scalar>("Mvap", 18.02)),
    Mmix_(dict.lookupOrDefault<scalar>("Mmix", 28.96))
{

    fvPatchField<scalar>::operator=
    (
        Field<scalar>("value", dict, p.size())
    );

    this->refValue() = *this;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = SMALL;
}

Foam::relativeHumidityInletOutletFvPatchScalarField::
relativeHumidityInletOutletFvPatchScalarField
(
    const relativeHumidityInletOutletFvPatchScalarField& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    relativeHumidity_(ptf.relativeHumidity_),
    phiName_(ptf.phiName_),
    TName_(ptf.TName_),
    PName_(ptf.PName_),
    rhoName_(ptf.rhoName_),
    Mvap_(ptf.Mvap_),
    Mmix_(ptf.Mmix_)
{}


Foam::relativeHumidityInletOutletFvPatchScalarField::
relativeHumidityInletOutletFvPatchScalarField
(
    const relativeHumidityInletOutletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    relativeHumidity_(ptf.relativeHumidity_),
    phiName_(ptf.phiName_),
    TName_(ptf.TName_),
    PName_(ptf.PName_),
    rhoName_(ptf.rhoName_),
    Mvap_(ptf.Mvap_),
    Mmix_(ptf.Mmix_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::relativeHumidityInletOutletFvPatchScalarField::patchDensity()
{
    const fvPatchScalarField& rho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    return tmp<scalarField>
    (
        new scalarField(rho)
    );

}

void Foam::relativeHumidityInletOutletFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchField<scalar>::autoMap(m);
    m(relativeHumidity_, relativeHumidity_);
}


void Foam::relativeHumidityInletOutletFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);

    const relativeHumidityInletOutletFvPatchScalarField& thftptf
        = refCast<const relativeHumidityInletOutletFvPatchScalarField>(ptf);
    relativeHumidity_.rmap(thftptf.relativeHumidity_, addr);
}


void Foam::relativeHumidityInletOutletFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    mixedFvPatchField<scalar>::autoMapGIB(mapper);
    mapper.map(relativeHumidity_, scalar(0));
}


void Foam::relativeHumidityInletOutletFvPatchScalarField::updateCoeffs()
{
    if (this->updated() || this->size() == 0)
    {
        return;
    }

    tmp<scalarField> rhop = patchDensity();

    // Get absolute pressure
    /**********************/
    scalarField pp(patch().lookupPatchField<volScalarField, scalar>(PName_));

    if (db().foundObject<dictionary>("transportProperties"))
    {
        // convert to [Pa]
        pp *= rhop();
        const IOdictionary& transportProperties =
            db().lookupObject<IOdictionary>("transportProperties");

        // pRef is already in [Pa]
        dimensionedScalar pRef("pRef",transportProperties.lookup("pRef"));
        pp += pRef.value();
    }
    else if (db().foundObject<dictionary>(basicThermo::dictName))
    {
        //check on compressible case
        //pRef = 0, do nothing
        const IOdictionary& materialProperties =
            db().lookupObject<IOdictionary>(basicThermo::dictName);

        // pRef is already in [Pa]
        scalar pRef = 0;
        if (materialProperties.isDict("referenceFields"))
        {
            const dictionary& refDict =
                materialProperties.subDict("referenceFields");
            pRef =
                refDict.found("p")
              ? refDict.lookup<dimensionedScalar>("p").value()
              : 0;
        }
        else if (materialProperties.found("pRef"))
        {
            pRef =
                materialProperties.lookup<scalar>("pRef");
        }
        pp += pRef;
    }
    else
    {
        FatalErrorInFunction
                << "No valid model for pRef"
                << exit(FatalError);
    }

    // calc patch partial pressure
    /*****************************/
    const fvPatchField<scalar>& Tp =
        patch().lookupPatchField<volScalarField, scalar> (TName_);

    //hardcoded saturation pressure for water vapour
    scalarField Psat
    (
        2337 * Foam::exp(6879*(1/293.15 - 1/Tp) - 5.031*Foam::log(Tp/293.15))
    );

    scalarField pi ( relativeHumidity_ * Psat );


    // calculate patch vapour concentration
    /************************************/
    this->refValue() = pi * Mvap_ / (pi * Mvap_ + (pp - pi)*Mmix_);

    const Field<scalar>& phip = this->patch().lookupPatchField
    (
        phiName_,
        reinterpret_cast<const surfaceScalarField*>(0),
        reinterpret_cast<const scalar*>(0)
    );

    this->valueFraction() = 1.0 - pos0(phip);

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::relativeHumidityInletOutletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
    writeEntryIfDifferent<word>(os, "p", "p", PName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<scalar>(os, "Mvap", 18.02, Mvap_);
    writeEntryIfDifferent<scalar>(os, "Mmix", 28.96, Mmix_);
    this->relativeHumidity_.writeEntry("relativeHumidity", os) ;
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::relativeHumidityInletOutletFvPatchScalarField::operator=
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
        relativeHumidityInletOutletFvPatchScalarField
    );
}

// ************************************************************************* //
