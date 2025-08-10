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
    (c) 2010-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels/derivedFvPatchFields/lumpedServerOutletTemperature/lumpedServerOutletTemperatureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

scalar lumpedServerOutletTemperatureFvPatchScalarField::calculateTeff
(
    const  scalar dt,
    const scalar tau1,
    const scalar Tamb,
    const scalar dTit
)
{
    return (tau1/(tau1 + dt))*Teff_
        + (dt/(tau1 + dt))*(Tamb + (scalar(1.0)-lambda_)*dTit);
}

scalar lumpedServerOutletTemperatureFvPatchScalarField::calculateDeltaT
(
    const  scalar dt,
    const scalar tau1,
    const scalar tau2,
    const scalar Tamb,
    const scalar dTit
)
{
    return dTit + (tau2/(tau1 + dt))*(Teff_ - Tamb - (scalar(1.0)-lambda_)*dTit);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lumpedServerOutletTemperatureFvPatchScalarField::
lumpedServerOutletTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    sourcePatchRegExp_(),
    patchIDs_(0),
    curTimeIndex_(-1),
    phiName_("phi"),
    verbose_(false),
    Ms_(),
    Cps_(),
    effectiveness_(),
    lambda_(),
    Ps_(),
    qs_(),
    Teff_(),
    deltaT_()
{}

lumpedServerOutletTemperatureFvPatchScalarField::
lumpedServerOutletTemperatureFvPatchScalarField
(
    const lumpedServerOutletTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    sourcePatchRegExp_(ptf.sourcePatchRegExp_),
    patchIDs_(ptf.patchIDs_),
    curTimeIndex_(-1),
    phiName_(ptf.phiName_),
    verbose_(ptf.verbose_),
    Ms_(ptf.Ms_),
    Cps_(ptf.Cps_),
    effectiveness_(ptf.effectiveness_),
    lambda_(ptf.lambda_),
    Ps_(ptf.Ps_, false),
    qs_(ptf.qs_),
    Teff_(ptf.Teff_),
    deltaT_(ptf.deltaT_)
{}

lumpedServerOutletTemperatureFvPatchScalarField::
lumpedServerOutletTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    sourcePatchRegExp_(dict.lookup("patches")),
    patchIDs_(p.patch().boundaryMesh().patchSet(sourcePatchRegExp_)),
    curTimeIndex_(-1),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    verbose_(dict.lookupOrDefault<bool>("verbose", false)),
    Ms_(readScalar(dict.lookup("mass"))),
    Cps_(readScalar(dict.lookup("Cps"))),
    effectiveness_(readScalar(dict.lookup("effectiveness"))),
    lambda_(dict.lookupOrDefault<scalar>("lambda", 0.5)),
    Ps_(),
    qs_(dict.lookupOrDefault<scalar>("massFlowRate", 0.015)),
    Teff_(dict.lookupOrDefault<scalar>("Teff", 300.0)),
    deltaT_(dict.lookupOrDefault<scalar>("deltaT", 0.0))
{

    if (returnReduce(size(), sumOp<label>()))
    {
        if (dict.found("power"))
        {
            Ps_ = Function1<scalar>::New("power", dict);
        }
        else
        {
            FatalIOErrorInFunction
            (dict)
            << "Please supply 'power' in [W] for the server BC"
            << exit(FatalIOError);
        }
    }

}

lumpedServerOutletTemperatureFvPatchScalarField::
lumpedServerOutletTemperatureFvPatchScalarField
(
    const lumpedServerOutletTemperatureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    sourcePatchRegExp_(ptf.sourcePatchRegExp_),
    patchIDs_(ptf.patchIDs_),
    curTimeIndex_(-1),
    phiName_(ptf.phiName_),
    verbose_(ptf.verbose_),
    Ms_(ptf.Ms_),
    Cps_(ptf.Cps_),
    effectiveness_(ptf.effectiveness_),
    lambda_(ptf.lambda_),
    Ps_(ptf.Ps_, false),
    qs_(ptf.qs_),
    Teff_(ptf.Teff_),
    deltaT_(ptf.deltaT_)
{
}

lumpedServerOutletTemperatureFvPatchScalarField::
lumpedServerOutletTemperatureFvPatchScalarField
(
    const lumpedServerOutletTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    sourcePatchRegExp_(ptf.sourcePatchRegExp_),
    patchIDs_(ptf.patchIDs_),
    curTimeIndex_(-1),
    phiName_(ptf.phiName_),
    verbose_(ptf.verbose_),
    Ms_(ptf.Ms_),
    Cps_(ptf.Cps_),
    effectiveness_(ptf.effectiveness_),
    lambda_(ptf.lambda_),
    Ps_(ptf.Ps_, false),
    qs_(ptf.qs_),
    Teff_(ptf.Teff_),
    deltaT_(ptf.deltaT_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void lumpedServerOutletTemperatureFvPatchScalarField::updateCoeffs()
{

    if (updated() || (curTimeIndex_ == this->db().time().timeIndex()))
    {
        return;
    }

    curTimeIndex_ = this->db().time().timeIndex();

    //calculate average Tin from neighbor patch
    scalar totalArea = 0;
    scalar meanTin = pTraits<scalar>::zero;

    const volScalarField& cvf
    (
        this->db().lookupObject<volScalarField>(this->internalField().name())
    );

    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        label patchi = iter.key();

        totalArea += gSum(this->patch().boundaryMesh()[patchi].magSf());

        meanTin += gSum
        (
            this->patch().boundaryMesh()[patchi].magSf()
            * cvf.boundaryField()[patchi]
        );
    }

    meanTin /= totalArea;

    const compressible::turbulenceModel& turbModel =
        db().lookupObject<compressible::turbulenceModel>
        (
            IOobject::groupName
            (
                compressible::turbulenceModel::propertiesName,
                internalField().group()
             )
         );

    const label patchi = patch().index();

    const basicThermo& thermo =
        db().lookupObject<basicThermo>(basicThermo::dictName);

    const scalarField& pp = thermo.p().boundaryField()[patchi];
    const scalarField& pT = thermo.T().boundaryField()[patchi];

    const scalarField Cpf(thermo.Cp(pp, pT, patchi));
    const scalar Cpa(gAverage(Cpf));

    const scalar dt(db().time().deltaTValue());

    if (curTimeIndex_ > 1)
    {
        const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>(phiName_);
        const fvsPatchField<scalar>& phip =
            patch().patchField<surfaceScalarField, scalar>(phi);

        const tmp<volScalarField> trho = turbModel.rho();
        const volScalarField& rho = trho;
        const scalarField& rhop = rho.boundaryField()[patchi];

        qs_ = mag(gSum(phip*rhop));
    }

    const scalar t = db().time().timeOutputValue();
    const scalar dTit(Ps_->value(t)/(qs_*Cpa));

    // time constants [s]
    scalar tau2;
    scalar tau1;

    if (curTimeIndex_ == 1)
    {
        //initial condition
        tau2 = scalar(0.0);
        tau1 = scalar(0.0);

        // initial Teff old = Tamb +(1-lambda)*dTit
        Teff_ = calculateTeff(dt, tau1,  meanTin, dTit);
    }
    else
    {
        tau2 = Ms_*Cps_/(qs_*Cpa);
        tau1 = tau2/(effectiveness_);
    }

    deltaT_ = calculateDeltaT(dt, tau1, tau2, meanTin, dTit);

    // update Teff and stored as old for next update
    scalar TeffNew_ =  calculateTeff(dt, tau1,  meanTin, dTit);
    Teff_ = TeffNew_ ;

    scalarField& Tp = *this;
    Tp = meanTin + deltaT_;

    fixedValueFvPatchField<scalar>::updateCoeffs();

    if (verbose_)
    {
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " Tout [K]: " << gAverage(*this)
            << " Teff [K]:" << Teff_
            << " Tamb [K]:" << meanTin
            << " deltaT [K]:" << deltaT_
            << " cooling flow rate [kg/s]:" << qs_
            << " tau1 [s]:" << tau1
            << " tau2 [s]:" << tau2
            << endl;
    }
}


void lumpedServerOutletTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("patches", sourcePatchRegExp_);

    os.writeEntry("phi", phiName_);
    os.writeEntry("verbose", verbose_);
    os.writeEntry("mass", Ms_);
    os.writeEntry("Cps", Cps_);
    os.writeEntry("effectiveness", effectiveness_);
    os.writeEntry("lambda", lambda_);

    Ps_->writeData(os);

    os.writeEntry("massFlowRate", qs_);
    os.writeEntry("Teff", Teff_);
    os.writeEntry("deltaT", deltaT_);


    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    lumpedServerOutletTemperatureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
