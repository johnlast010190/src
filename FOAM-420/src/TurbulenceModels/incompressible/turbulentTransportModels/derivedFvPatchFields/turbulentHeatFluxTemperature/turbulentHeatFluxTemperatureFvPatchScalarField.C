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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "turbulentHeatFluxTemperatureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // declare specialization within 'Foam' namespace
    template<>
    const char* NamedEnum
    <
        Foam::incompressible::
        turbulentHeatFluxTemperatureFvPatchScalarField::heatSourceType,
        4
    >::names[] =
    {
        "power",
        "flux",
        "powerTV",
        "fluxTV"
    };

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const NamedEnum
<
    turbulentHeatFluxTemperatureFvPatchScalarField::heatSourceType,
    4
> turbulentHeatFluxTemperatureFvPatchScalarField::heatSourceTypeNames_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    heatSource_(hsPower),
    q_(p.size(), 0.0),
    qtimevar_()
{}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const turbulentHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    heatSource_(ptf.heatSource_),
    q_(mapper(ptf.q_)),
    qtimevar_(ptf.qtimevar_, false)
{}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    heatSource_(heatSourceTypeNames_.read(dict.lookup("heatSource"))),
    q_(p.size(), 0.0),
    qtimevar_()
{
    if (returnReduce(size(), sumOp<label>()))
    {
        if (dict.found("q"))
        {
            q_ = scalarField("q", dict, p.size());
        }
        else if (dict.found("qtime"))
        {
            qtimevar_ = Function1<scalar>::New("qtime", dict);
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Please supply either 'q' or"
                << "qtime" << exit(FatalIOError);
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
}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const turbulentHeatFluxTemperatureFvPatchScalarField& thftpsf
)
:
    fixedGradientFvPatchScalarField(thftpsf),
    heatSource_(thftpsf.heatSource_),
    q_(thftpsf.q_),
    qtimevar_(thftpsf.qtimevar_, false)
{}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const turbulentHeatFluxTemperatureFvPatchScalarField& thftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(thftpsf, iF),
    heatSource_(thftpsf.heatSource_),
    q_(thftpsf.q_),
    qtimevar_(thftpsf.qtimevar_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    m(q_, q_);
}


void turbulentHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const turbulentHeatFluxTemperatureFvPatchScalarField& thftptf =
        refCast<const turbulentHeatFluxTemperatureFvPatchScalarField>
        (
            ptf
        );

    q_.rmap(thftptf.q_, addr);
}


void turbulentHeatFluxTemperatureFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedGradientFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(q_, scalar(0));
}


void turbulentHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const incompressibleTurbulenceModel& turbModel =
        db().lookupObject<incompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    const tmp<volScalarField> talphaEff
        = turbModel.alphaEff()*turbModel.rho(); //<-- this rho is only used for incompressible
    const volScalarField& alphaEff( talphaEff );
    const scalarField& alphaEffp( alphaEff.boundaryField()[patchi] );

    const tmp<volScalarField> tCp( turbModel.Cp() );
    const volScalarField& Cp( tCp );
    const scalarField& Cpp( Cp.boundaryField()[patchi] );

    scalarField Qt(patch().size());

    switch (heatSource_)
    {
        case hsPower:
        {
            Qt = q_;
            const scalar Ap = gSum(patch().magSf());
            Qt /= Ap;
            break;
        }
        case hsFlux:
        {
            Qt = q_;
            break;
        }
        case hsTVPower:
        {
            const scalar t = db().time().timeOutputValue();
            scalar qtimevar = qtimevar_->value(t);
            Qt = qtimevar;

            const scalar Ap = gSum(patch().magSf());
                Qt /= Ap;
            break;
        }
        case hsTVFlux:
        {
            const scalar t = db().time().timeOutputValue();
            scalar qtimevar = qtimevar_->value(t);
            Qt = qtimevar;
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown heat source type. Valid types are: "
                << heatSourceTypeNames_ << nl << exit(FatalError);
        }
    }


    // Update option sources, using Newton's method to account for
    // possible face value dependence
    // We are solving
    // kappaEff*deltaCoeff*(Twall-Tpi) = Qt + boundarySource(Twall)
    // or f := Qt + boundarySource(Twall) - kappaEff*deltaCoeffs*(Twall-Tpi) = 0

    scalarField Tpi( this->patchInternalField() );
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    scalarField C1( (Qt + Cpp*alphaEffp*Tpi*deltaCoeffs) );
    scalarField C2( Cpp*alphaEffp*deltaCoeffs );

    label i = 0;
    scalar Tw_err = GREAT;
    scalar maxErr = 1e-5;
    label maxLoops = 100;

    scalarField Tw_new = *this;
    scalarField Tw_old;

    scalarField Qopt;
    do
    {
        Tw_old = Tw_new;
        scalarField bSourceDeriv;
        Qopt = boundarySources(Tw_old, bSourceDeriv);

        scalarField f( C1 + Qopt - C2*Tw_old );
        scalarField df( bSourceDeriv - C2 );

        Tw_new = Tw_old - f/df;
        Tw_err = gMax(mag(Tw_new-Tw_old)/stabilise(Tw_old, SMALL));
        i++;
    }
    while (i < maxLoops && Tw_err > maxErr);

    if (i == maxLoops)
    {
        WarningInFunction
            << "Non-convergence in Newton's method for "
            << "temperature-dependent boundary source, patch: "
            << this->patch().name() << nl << endl;
    }

    gradient() = (Qt+Qopt)/(Cpp*alphaEffp);
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void turbulentHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("heatSource", heatSourceTypeNames_[heatSource_]);
    gradient().writeEntry("gradient", os);

    switch (heatSource_)
    {
        case hsPower:
        {
            q_.writeEntry("q", os);
            break;
        }
        case hsFlux:
        {
            q_.writeEntry("q", os);
            break;
        }
        case hsTVPower:
        {
            qtimevar_->writeData(os);
            break;
        }
        case hsTVFlux:
        {
            qtimevar_->writeData(os);
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown heat source type. Valid types are: "
                << heatSourceTypeNames_ << nl << exit(FatalError);
        }
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentHeatFluxTemperatureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
