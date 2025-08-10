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
    (c) 1991-2009 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "turbulentMassFluxConcentrationFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "turbulenceModel.H"
#include "basicThermo/basicThermo.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // declare specialization within 'Foam' namespace
    template<>
    const char* NamedEnum
    <
        turbulentMassFluxConcentrationFvPatchScalarField::massFluxType,
        2
    >::names[] =
    {
        "massflux",
        "timemassflux"
    };

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const NamedEnum
<
   turbulentMassFluxConcentrationFvPatchScalarField::massFluxType,
   2
>  turbulentMassFluxConcentrationFvPatchScalarField::massFluxTypeNames_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentMassFluxConcentrationFvPatchScalarField::
turbulentMassFluxConcentrationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    m_(p.size(), 0.0),
    mt_(),
    speciesName_("dummy"),
    rhoName_("rho"),
    Dlam_(dimensionedScalar("D", sqr(dimLength)/dimTime, 0.0)),
    massflux_(massFluxS)
{}


turbulentMassFluxConcentrationFvPatchScalarField::
turbulentMassFluxConcentrationFvPatchScalarField
(
    const turbulentMassFluxConcentrationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    m_(mapper(ptf.m_)),
    mt_(ptf.mt_, false),
    speciesName_(ptf.speciesName_),
    rhoName_(ptf.rhoName_),
    Dlam_(ptf.Dlam_),
    massflux_(ptf.massflux_)
{}


turbulentMassFluxConcentrationFvPatchScalarField::
turbulentMassFluxConcentrationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    m_(p.size(), 0.0),
    mt_(),
    speciesName_(fvPatchScalarField::internalField().name()),
    rhoName_(dict.lookupOrDefault<word>("rho","rho")),
    Dlam_(dimensionedScalar("D"+speciesName_, sqr(dimLength)/dimTime, 0.0)),
    massflux_(massFluxTypeNames_.read(dict.lookup("massFluxType")))
{
    //strip trailing "_0" from old time field species name
    string underscore0("_0");
    string sn(speciesName_);

    if (sn.find(underscore0, sn.size()-2) != string::npos)
    {
        //delete last two chars
        speciesName_ = sn.substr(0, sn.size() - underscore0.size());
    }

    if (db().foundObject<IOdictionary>("transportProperties"))
    {
        const IOdictionary& transportProperties
            = db().lookupObject<IOdictionary>("transportProperties");

        Dlam_ = dimensionedScalar
        (
            "D"+speciesName_,
            transportProperties.lookup("D"+speciesName_)
        );
    }
    else if (db().foundObject<dictionary>(basicThermo::dictName))
    {
        const IOdictionary& thermophysicalProperties =
            db().lookupObject<IOdictionary>(basicThermo::dictName);

        if (basicThermo::dictName == basicThermo::matDictName)
        {
            const word subDictName(speciesName_ + "concentrationTransport");
            Dlam_ = dimensionedScalar
            (
                "D",
                thermophysicalProperties.optionalSubDict(subDictName).lookup("D")
            );

        }
        else
        {
            Dlam_ = dimensionedScalar
            (
                "D" + speciesName_,
                thermophysicalProperties.lookup("D" + speciesName_)
            );
        }
    }
    else
    {
        WarningInFunction
            << "Laminar diffusivity for " << speciesName_
            << " not initialised. This message can be ignored during "
            << "initialisation." << endl;
    }

    if (dict.found("massflux"))
    {
        m_ = scalarField("massflux", dict, p.size());
    }
    else if (dict.found("timemassflux"))
    {
            mt_ = Function1<scalar>::New("timemassflux", dict);
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Please supply either 'massflux' or"
            << "timemassflux" << exit(FatalIOError);
    }

    gradient() = 0.0;

    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        }

    Field<scalar>::operator=
    (
        this->patchInternalField() + gradient()/this->patch().deltaCoeffs()
    );
}


turbulentMassFluxConcentrationFvPatchScalarField::
turbulentMassFluxConcentrationFvPatchScalarField
(
    const turbulentMassFluxConcentrationFvPatchScalarField& thftpsf
)
:
    fixedGradientFvPatchScalarField(thftpsf),
    m_(thftpsf.m_),
    mt_(thftpsf.mt_, false),
    speciesName_(thftpsf.speciesName_),
    Dlam_(thftpsf.Dlam_),
    massflux_(thftpsf.massflux_)
{}


turbulentMassFluxConcentrationFvPatchScalarField::
turbulentMassFluxConcentrationFvPatchScalarField
(
    const turbulentMassFluxConcentrationFvPatchScalarField& thftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(thftpsf, iF),
    m_(thftpsf.m_),
    mt_(thftpsf.mt_, false),
    speciesName_(thftpsf.speciesName_),
    rhoName_(thftpsf.rhoName_),
    Dlam_(thftpsf.Dlam_),
    massflux_(thftpsf.massflux_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentMassFluxConcentrationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    m(m_, m_);
}


void turbulentMassFluxConcentrationFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const turbulentMassFluxConcentrationFvPatchScalarField& thftptf =
        refCast<const turbulentMassFluxConcentrationFvPatchScalarField>
        (
            ptf
        );

    m_.rmap(thftptf.m_, addr);
}


void turbulentMassFluxConcentrationFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedGradientFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(m_, scalar(0));
}


void turbulentMassFluxConcentrationFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchI = patch().index();

    const volScalarField& Dt
        = db().lookupObject<volScalarField>("Dt"+speciesName_);

    scalarField Deffp ( Dlam_.value() + scalarField(Dt.boundaryField()[patchI]) );

    //density on patch
    scalarField rhoap(this->size(),0.0);

    if (db().foundObject<volScalarField>(rhoName_))
    {
        rhoap =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);
    }
    else if (db().foundObject<dictionary>("transportProperties"))
    {
        const IOdictionary& transportProperties
            = db().lookupObject<IOdictionary>("transportProperties");

        dimensionedScalar rhoRef("rhoRef",transportProperties.lookup("rho"));
        rhoap = rhoRef.value();
    }
    else
    {
        FatalErrorInFunction
            << "No access to rho"
            << exit(FatalError);
    }

    switch (massflux_)
    {
        case massFluxS:
    {
        gradient() = m_/(rhoap*Deffp);
        break;
    }
        case massFluxTV:
    {
        const scalar t = db().time().timeOutputValue();
        gradient() = (mt_->value(t))/(rhoap*Deffp);
        break;
    }
        default:
        {
            FatalErrorInFunction
                << "Unknown mass flux type. Valid types are: "
                << massFluxTypeNames_ << nl << exit(FatalError);
        }
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void turbulentMassFluxConcentrationFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("massFluxType", massFluxTypeNames_[massflux_]);

    switch (massflux_)
    {
        case massFluxS:
    {
        m_.writeEntry("massflux", os);
        break;
    }
        case massFluxTV:
    {
        mt_->writeData(os);
        break;
    }
        default:
        {
            FatalErrorInFunction
                << "In write, unknown mass flux type. Valid types are: "
                << massFluxTypeNames_ << nl << exit(FatalError);
        }
    }

    gradient().writeEntry("gradient", os);
    writeEntry("value", os);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentMassFluxConcentrationFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
