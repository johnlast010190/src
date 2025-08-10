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
    (c) 1991-2008 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "phaseChangeHumidityFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "turbulenceModel.H"
#include "RAS/RASModel/RASModel.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/ddtSchemes/steadyStateDdtScheme/steadyStateDdtScheme.H"
#include "basicThermo/basicThermo.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * *  Private member functions * * * * * * * * * * * //

bool phaseChangeHumidityFvPatchScalarField::timeScheme()
{
    word ddtScheme
    (
        this->internalField().mesh()
       .schemes().ddtScheme(this->internalField().name())
    );

    if
    (
        ddtScheme == fv::steadyStateDdtScheme<scalar>::typeName
    )
    {
        return false;
    }
    else
    {
        return true;
    }
}

tmp<scalarField> phaseChangeHumidityFvPatchScalarField::massDiffusionCoeff()
{
    //read laminar mass diffusivity
    dimensionedScalar D("D", dimArea/dimTime, 0.0);

    if (db().foundObject<dictionary>("transportProperties"))
    {
        const IOdictionary& transportProperties
            = db().lookupObject<IOdictionary>("transportProperties");

        D = dimensionedScalar("Dw",transportProperties.lookup("Dw"));
    }
    else if (db().foundObject<dictionary>(basicThermo::dictName))
    {
        const IOdictionary& thermophysicalProperties
            = db().lookupObject<IOdictionary>(basicThermo::dictName);

        bool isMat = (basicThermo::dictName == basicThermo::matDictName);
        // In general case it isn't clear what the sub-dict should be called
        if (isMat)
        {
            //TODO Should this be a name of the sub-dict?
            D =
                dimensionedScalar
                (
                    "D",
                    thermophysicalProperties.optionalSubDict
                    (
                        "wconcentrationTransport"
                    ).lookup("D")
                );
        }
        else
        {
            D = dimensionedScalar("Dw",thermophysicalProperties.lookup("Dw"));
        }
    }
    else
    {
        FatalErrorInFunction
            << "No access to Dw"
            << exit(FatalError);
    }

    //lookup turbulent mass diffusivity
    const fvPatchField<scalar>& Dtp =
        patch().lookupPatchField<volScalarField, scalar>(DtName_);

    return tmp<scalarField>
    (
        new scalarField
        (
            (Dtp + D.value())*wetRatio_
        )
    );
}

tmp<scalarField> phaseChangeHumidityFvPatchScalarField::patchDensity()
{
    const fvPatchScalarField& rho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    return tmp<scalarField>
    (
        new scalarField(rho)
    );
}

scalar phaseChangeHumidityFvPatchScalarField::dropDiameter(label i, scalar theta)
{
    scalar H = filmMass_[i] / rhoLiquid_;

    if (theta > SMALL)
    {
        return
        (
            12.0*sqrt(3.0)*sqr(sin(theta))*H
            / constant::mathematical::pi
            / sqr(1.0 - cos(theta))
            /(2.0 + cos(theta))
        );
    }
    else
    {
        return SMALL;
    }
}

scalar phaseChangeHumidityFvPatchScalarField::contactAngle(label i, scalar d)
{
    scalar H = filmMass_[i] / rhoLiquid_;

    scalar theta = thetaA_;

    //Newton-Raphson solver for theta
    label j = 0;
    scalar theta_err = GREAT;
    label maxLoops = 100;
    scalar maxErr = 0.01;

    do
    {

        scalar C1 = 12.0*sqrt(3.0)*H / constant::mathematical::pi;
        scalar sint = sin(theta);
        scalar cost = cos(theta);
        scalar f1 = sqr(sint);
        scalar f2 = 1/sqr(1.0 - cost);
        scalar f3 = (2.0 + cost);

        scalar f = C1 * f1 * f2 - d * f3;

        scalar df
            = C1*(2.0*sint*cost*f2 - 2.0*sint/pow3(1.0-cost)*f1)
            - sint*d;


        scalar theta_old = theta;

        theta = theta_old - f/(sign(df)*max(mag(df), SMALL));
        theta_err = 2*mag(theta-theta_old)/max(SMALL, (theta+theta_old));
        j++;

        if (theta < SMALL)
        {
            theta = SMALL;
            break;
        }

    } while (j < maxLoops && theta_err > maxErr);

    if (j==maxLoops)
    {
        Info<< j << " exceeded max loops for theta calculation. theta : " << theta
            << ", d: " << d << ", H : " << H <<
             endl;

    }

    if (theta > (thetaA_+SMALL))
    {
        WarningInFunction
            << "Error in contact angle prediction."
            << "theta : " << theta << " > advancing contact angle : "
            << thetaA_ << endl;
    }

    theta = max(theta, SMALL);
    return theta;
}

void phaseChangeHumidityFvPatchScalarField::updateWetDryAreaRatio()
{
    //wetRatio_ =  dryWetCondensing_;

    scalarField wpi ( this->patchInternalField() );
    const scalarField& wp = *this;
    scalar pie = constant::mathematical::pi;


    forAll(wp, i)
    {

        if (filmMass_[i] < VSMALL)
        {
            wetRatio_[i] = 1;
        }
        else if (wp[i] > wpi[i]) //evaporation
        {

            scalar ddrop = dDropCondensing_()[i];
            scalar theta = contactAngle(i, ddrop);

            scalar ndrops
                = dryWetCondensing_*patch().magSf()[i]/(0.25*pie*sqr(ddrop));

            if (theta < thetaR_)
            {
                theta = thetaR_;
                ddrop = dropDiameter(i, theta);
            }

            scalar Acap
                 = ndrops * 0.5 * pie * sqr(ddrop)
                 * (1 - cos(theta))/sqr(sin(theta));

            wetRatio_[i] = Acap / patch().magSf()[i];

        }
        else // condensing
        {
            scalar theta = thetaA_;
            scalar ddrop = dropDiameter(i, theta);
            dDropCondensing_()[i] = ddrop;

            scalar ndrops
                = dryWetCondensing_*patch().magSf()[i]/(0.25*pie*sqr(ddrop));

            scalar Acap
                 = ndrops * 0.5 * pie * sqr(ddrop)
                 * (1 - cos(theta))/sqr(sin(theta));

            wetRatio_[i] = Acap / patch().magSf()[i];

        }
    }
}





// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

phaseChangeHumidityFvPatchScalarField::phaseChangeHumidityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    TName_("T"),
    PName_("p"),
    rhoName_("rho"),
    DtName_("Dtw"),
    filmMass_(p.size()),
    heatFlux_(p.size()),
    wetRatio_(p.size(), 1.0),
    transient_(false),
    steadyEvaporation_(false),
    thermalCoupling_(false),
    dHv_(2.257e+6),
    Mvap_(18.02),
    Mmix_(28.96),
    dropSim_(false),
    dDropCondensing_(nullptr),
    thetaA_(90.0),
    thetaR_(10.0),
    dryWetCondensing_(0.55),
    rhoLiquid_(1000.0)

{
    fvPatchScalarField::operator= (pTraits<scalar>::zero);
    this->refValue() = pTraits<scalar>::zero ;
    this->refGrad() = pTraits<scalar>::zero ;
    this->valueFraction() = 0.0 ;
}

phaseChangeHumidityFvPatchScalarField::phaseChangeHumidityFvPatchScalarField
(
    const phaseChangeHumidityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    PName_(ptf.PName_),
    rhoName_(ptf.rhoName_),
    DtName_(ptf.DtName_),
    filmMass_(mapper(ptf.filmMass_)),
    heatFlux_(mapper(ptf.heatFlux_)),
    wetRatio_(mapper(ptf.wetRatio_)),
    transient_(ptf.transient_),
    steadyEvaporation_(ptf.steadyEvaporation_),
    thermalCoupling_(ptf.thermalCoupling_),
    dHv_(ptf.dHv_),
    Mvap_(ptf.Mvap_),
    Mmix_(ptf.Mmix_),
    dropSim_(ptf.dropSim_),
    dDropCondensing_(),
    thetaA_(ptf.thetaA_),
    thetaR_(ptf.thetaR_),
    dryWetCondensing_(ptf.dryWetCondensing_),
    rhoLiquid_(ptf.rhoLiquid_)
{
    if (ptf.dDropCondensing_.valid())
    {
        dDropCondensing_.reset(mapper(ptf.dDropCondensing_()).ptr());
    }
}

phaseChangeHumidityFvPatchScalarField::phaseChangeHumidityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF, dict, false),
    TName_(dict.lookupOrDefault<word>("T","T")),
    PName_(dict.lookupOrDefault<word>("p","p")),
    rhoName_(dict.lookupOrDefault<word>("rho","rho")),
    DtName_(dict.lookupOrDefault<word>("Dtw","Dtw")),
    filmMass_("filmMass", dict, p.size()),
    heatFlux_(p.size(), 0.0),
    wetRatio_(p.size(), 1.0),
    transient_(false),
    steadyEvaporation_
    (
        dict.lookupOrDefault<Switch>("steadyEvaporation", false)
    ),
    thermalCoupling_
    (
        dict.lookupOrDefault<Switch>("thermalCoupling", false)
    ),
    dHv_(dict.lookupOrDefault<scalar>("heatOfEvaporation", 2.257e+6)),
    Mvap_(dict.lookupOrDefault<scalar>("Mvap", 18.02)),
    Mmix_(dict.lookupOrDefault<scalar>("Mmix", 28.96)),
    dropSim_(dict.lookupOrDefault<Switch>("droplets", false)),
    dDropCondensing_(nullptr),
    thetaA_(-1),
    thetaR_(-1),
    dryWetCondensing_(-1),
    rhoLiquid_(-1)
{
    this->refGrad() = 0;
    this->valueFraction() = 0.0 ;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
        refValue() = *this;
        this->valueFraction() = 1.0 ;
    }
    else
    {
        mixedFvPatchField<scalar>::updateCoeffs();
    }

    if (dict.found("heatFlux"))
    {
        heatFlux_ = scalarField("heatFlux", dict, p.size());
    }

    if (dropSim_)
    {
        scalar pi180 = constant::mathematical::pi/180.0;
        thetaA_ = pi180*readScalar(dict.lookup("contactAngle"));
        thetaR_ = pi180*dict.lookupOrDefault<scalar>("retreatingContactAngle", 10.0);
        dryWetCondensing_ = dict.lookupOrDefault<scalar>("liquidAreaRatio", 0.55);
        rhoLiquid_ = dict.lookupOrDefault<scalar>("rhoLiquid", 1000.0);

        // set initial droplet size for evaporating interfaces
        if (dict.found("condensingDropletDiameter"))
        {
            dDropCondensing_.reset
            (
                new scalarField("condensingDropletDiameter", dict, p.size())
            );
        }
        else
        {
            dDropCondensing_.set(new scalarField(p.size(), 0));

            forAll(dDropCondensing_(), di)
            {
                dDropCondensing_()[di] = dropDiameter(di, thetaA_);
            }
        }
        updateWetDryAreaRatio();
    }
}

phaseChangeHumidityFvPatchScalarField::phaseChangeHumidityFvPatchScalarField
(
  const phaseChangeHumidityFvPatchScalarField& ptf
)
    :
    mixedFvPatchField<scalar>(ptf),
    TName_(ptf.TName_),
    PName_(ptf.PName_),
    rhoName_(ptf.rhoName_),
    DtName_(ptf.DtName_),
    filmMass_(ptf.filmMass_),
    heatFlux_(ptf.heatFlux_),
    wetRatio_(ptf.wetRatio_),
    transient_(ptf.transient_),
    steadyEvaporation_(ptf.steadyEvaporation_),
    thermalCoupling_(ptf.thermalCoupling_),
    dHv_(ptf.dHv_),
    Mvap_(ptf.Mvap_),
    Mmix_(ptf.Mmix_),
    dropSim_(ptf.dropSim_),
    dDropCondensing_(),
    thetaA_(ptf.thetaA_),
    thetaR_(ptf.thetaR_),
    dryWetCondensing_(ptf.dryWetCondensing_),
    rhoLiquid_(ptf.rhoLiquid_)
{
    if (ptf.dDropCondensing_.valid())
    {
        dDropCondensing_.reset(new scalarField(ptf.dDropCondensing_()));
    }
}

phaseChangeHumidityFvPatchScalarField::phaseChangeHumidityFvPatchScalarField
(
  const phaseChangeHumidityFvPatchScalarField& ptf,
  const DimensionedField<scalar, volMesh>& iF
)
    :
    mixedFvPatchField<scalar>(ptf, iF),
    TName_(ptf.TName_),
    PName_(ptf.PName_),
    rhoName_(ptf.rhoName_),
    DtName_(ptf.DtName_),
    filmMass_(ptf.filmMass_),
    heatFlux_(ptf.heatFlux_),
    wetRatio_(ptf.wetRatio_),
    transient_(ptf.transient_),
    steadyEvaporation_(ptf.steadyEvaporation_),
    thermalCoupling_(ptf.thermalCoupling_),
    dHv_(ptf.dHv_),
    Mvap_(ptf.Mvap_),
    Mmix_(ptf.Mmix_),
    dropSim_(ptf.dropSim_),
    dDropCondensing_(),
    thetaA_(ptf.thetaA_),
    thetaR_(ptf.thetaR_),
    dryWetCondensing_(ptf.dryWetCondensing_),
    rhoLiquid_(ptf.rhoLiquid_)
{
    if (ptf.dDropCondensing_.valid())
    {
        dDropCondensing_.reset(new scalarField(ptf.dDropCondensing_()));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void phaseChangeHumidityFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    m(filmMass_, filmMass_);
    m(heatFlux_, heatFlux_);
    m(wetRatio_, wetRatio_);

    if (dDropCondensing_.valid())
    {
        m(*dDropCondensing_, *dDropCondensing_);
    }
}


void phaseChangeHumidityFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const phaseChangeHumidityFvPatchScalarField& thftptf
        = refCast<const phaseChangeHumidityFvPatchScalarField>(ptf);

    filmMass_.rmap(thftptf.filmMass_, addr);
    heatFlux_.rmap(thftptf.heatFlux_, addr);
    wetRatio_.rmap(thftptf.wetRatio_, addr);

    if (dDropCondensing_.valid())
    {
        dDropCondensing_->rmap(thftptf.dDropCondensing_(), addr);
    }
}


void phaseChangeHumidityFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(filmMass_, scalar(0));
    mapper.map(heatFlux_, scalar(0));
    mapper.map(wetRatio_, scalar(0));
    if (dDropCondensing_.valid())
    {
        mapper.map(dDropCondensing_(), scalar(0));
    }
}


void phaseChangeHumidityFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    tmp<scalarField> rhop = patchDensity();

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

    // calc wall saturation pressure
    /*******************************/
    const fvPatchField<scalar>& Tp =
        patch().lookupPatchField<volScalarField, scalar> (TName_);

    //hardcoded saturation pressure for water vapour
    scalarField Psat
    (
        2337 * Foam::exp(6879*(1/293.15 - 1/Tp) - 5.031*Foam::log(Tp/293.15))
    );


    // calculate wall vapour concentration
    /************************************/
    refValue() = Psat * Mvap_ / (Psat * Mvap_ + (pp - Psat)*Mmix_);
    valueFraction() = 1.0;

    // if evaporating, check if water present, if not set to zero gradient
    /********************************************************************/
    scalarField wc ( this->patchInternalField() );

    transient_ = timeScheme();

    forAll(wc, i)
    {
        if (transient_)
        {
            if (refValue()[i] >= wc[i] && filmMass_[i] < VSMALL)
            {
                // wants to evaporate, but there is no liquid on wall
                valueFraction()[i] = 0.0;
            }
        }
        else
        {
            if (refValue()[i] >= wc[i] && !steadyEvaporation_)
            {
                // wants to evaporate, but no steady evaporation allowed
                valueFraction()[i] = 0.0;
            }
        }
    }

    // need to update boundary value here to calculate mass and heat fluxes
    fvPatchScalarField::operator=
    (
        valueFraction()*refValue()
      +
        (1.0 - valueFraction())*
        (
            this->patchInternalField()
          + refGrad() /this->patch().deltaCoeffs()
        )
    );
    mixedFvPatchField<scalar>::updateCoeffs();


    /***********************************************************************/
    tmp<scalarField> Deff = massDiffusionCoeff();

    scalarField massFlowRate ( - rhop() * Deff() * this->snGrad() );

    if (transient_)
    {
        filmMass_ +=  db().time().deltaT().value()* massFlowRate;

        // limit minimum mass to zero
        filmMass_ = max(filmMass_, scalar(0.0));
    }
    else
    {
        filmMass_ = massFlowRate;
    }

    if (thermalCoupling_)
    {
        // to get rid of the relaxation, the new temperature has to be
        // made more implicit
        // probably a good idea for transient runs
        scalar relax = 0.1;
        heatFlux_ = heatFlux_ * (1- relax) + relax * massFlowRate * dHv_;
    }
    else
    {
        heatFlux_ = 0;
    }

    if (dropSim_)
    {
        updateWetDryAreaRatio();
    }

    /*
    Info<< "heatFlux_ : " << gAverage(heatFlux_)
         << ", massFlowRate : " << gAverage(massFlowRate) << endl;
    Info<< "t ds :" << transient_ << dropSim_ << endl;
    */

    const fvMesh& mesh = this->internalField().mesh();
    const label patchi = patch().index();

    // store and write filmMass field
    if (db().found("filmMass"))
    {
        const volScalarField& fmc =
            db().lookupObject<volScalarField>("filmMass");

        volScalarField& fm = const_cast<volScalarField&>(fmc);

        fm.boundaryFieldRef()[patchi] = filmMass();
    }
    else
    {
        db().store
        (
            new volScalarField
            (
                IOobject
                (
                    "filmMass",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "filmMass",
                    dimMass/dimArea,
                    0.0
                )
            )
        );

        const volScalarField& fmc =
            db().lookupObject<volScalarField>("filmMass");

        volScalarField& fm = const_cast<volScalarField&>(fmc);

        fm.boundaryFieldRef()[patchi] = filmMass();
    }


}

void phaseChangeHumidityFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os) ;
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
    writeEntryIfDifferent<word>(os, "p", "P", PName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "Dtw", "Dtw", DtName_);
    writeEntryIfDifferent<scalar>(os, "heatOfEvaporation",  2.257e+6, dHv_);
    writeEntryIfDifferent<scalar>(os, "Mvap", 18.02, Mvap_);
    writeEntryIfDifferent<scalar>(os, "Mmix", 28.96, Mmix_);
    writeEntryIfDifferent<Switch>
    (
        os, "steadyEvaporation", false, steadyEvaporation_
    );
    writeEntryIfDifferent<Switch>
    (
        os, "thermalCoupling", false, thermalCoupling_
    );

    writeEntryIfDifferent<Switch>(os, "droplets", false, dropSim_);

    if (dropSim_)
    {
        scalar pi180 = constant::mathematical::pi/180.0;

        os.writeEntry("contactAngle", (thetaA_/pi180));

        writeEntryIfDifferent<scalar>
        (
            os, "retreatingContactAngle", 10.0, thetaR_/pi180
        );
        writeEntryIfDifferent<scalar>
        (
            os, "liquidAreaRatio", 0.55, dryWetCondensing_
        );
        writeEntryIfDifferent<scalar>
        (
            os, "rhoLiquid", 1000.0, rhoLiquid_
        );
        if (dDropCondensing_.valid())
        {
            dDropCondensing_->writeEntry("condensingDropletDiameter", os);
        }
    }
    filmMass_.writeEntry("filmMass", os) ;
    heatFlux_.writeEntry("heatFlux", os);
    this->writeEntry("value", os) ;
}

makePatchTypeField(fvPatchScalarField, phaseChangeHumidityFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
