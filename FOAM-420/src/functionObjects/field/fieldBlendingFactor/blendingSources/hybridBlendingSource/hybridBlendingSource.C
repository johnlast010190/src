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
    (c) 2011 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "fieldBlendingFactor/blendingSources/hybridBlendingSource/hybridBlendingSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "basicThermo/basicThermo.H"
#include "fields/fvsPatchFields/basic/fixedValue/fixedValueFvsPatchFields.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"

// Use surfaceInterpolate from finiteVolume, not this library!
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "finiteVolume/fvc/fvcAverage.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace blendingSources
{

defineTypeNameAndDebug(hybridBlendingSource, 0);
addToRunTimeSelectionTable
(blendingSource, hybridBlendingSource, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- return effective kinematic viscosity
tmp<volScalarField> hybridBlendingSource::nuEff()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            obr_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        if (isA<compressible::LESModel>(turb))
        {
            return (turb.muEff()/turb.rho());
        }
    }
    else if (mesh_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const icoTurbModel& turb =
            obr_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        if (isA<incompressible::LESModel>(turb))
        {
            return (turb.nuEff());
        }
    }

    FatalError << this->typeName << " blendingSource"
               << " requires an LES model to function."
               << exit(FatalError);

    return tmp<volScalarField>();

}

//- return velocity
tmp<volVectorField> hybridBlendingSource::U()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;


    if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            obr_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        if (isA<compressible::LESModel>(turb))
        {
            return tmp<volVectorField>(turb.U());
        }
    }
    else if (mesh_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const icoTurbModel& turb =
            obr_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        if (isA<incompressible::LESModel>(turb))
        {
            return tmp<volVectorField>(turb.U());
        }
    }

    FatalError << this->typeName << " blendingSource"
               << " requires an LES model to function."
               << exit(FatalError);

    return tmp<volVectorField>();

}

//- return CDES*delta
tmp<volScalarField> hybridBlendingSource::delta()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            obr_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        if (isA<compressible::LESModel>(turb))
        {
            return refCast<const compressible::LESModel>(turb).delta();
        //    return (tmp<volScalarField>(turb.delta()));
        }
    }
    else if (mesh_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const icoTurbModel& turb =
            obr_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        if (isA<incompressible::LESModel>(turb))
        {
            return refCast<const incompressible::LESModel>(turb).delta();
            //return (tmp<volScalarField>(turb.delta()));
        }
    }

    FatalError << this->typeName << " blendingSource"
               << " requires an LES model to function."
               << exit(FatalError);

    return tmp<volScalarField>();

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hybridBlendingSource::hybridBlendingSource
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    blendingSource(obr, mesh, dict),
    lRef_(dict.lookupOrDefault<scalar>("lRef", -1)),
    Uref_(dict.lookupOrDefault<scalar>("Uref", -1)),
    CH1_(dict.lookupOrDefault<scalar>("CH1", 3)),
    CH2_(dict.lookupOrDefault<scalar>("CH2", 1)),
    CH3_(dict.lookupOrDefault<scalar>("CH3", 2)),
    CDES_(dict.lookupOrDefault<scalar>("CDES", 0.65)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    OmegaLim_(dict.lookupOrDefault<scalar>("omegaLim", 1e-03)),
    sigmaMax_(dict.lookupOrDefault<scalar>("sigmaMax", 1)),
    sigmaMin_(dict.lookupOrDefault<scalar>("sigmaMin", 0)),
    approxTanh_(dict.lookupOrDefault<Switch>("approxTanh", true))
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField> hybridBlendingSource::sourceField()
{
    tmp<volTensorField> gradU(fvc::grad(U()));

    tmp<volScalarField> Ostar
    (
        sqrt(2.0)*mag(skew(gradU()))
    );

    tmp<volScalarField> Sstar
    (
        sqrt(2.0)*mag(symm(gradU()))
    );
    gradU.clear();

    if (lRef_ == -1)
        lRef_ = cbrt(gSum(mesh_.V()));
    if (Uref_ == -1)
        Uref_ = max(SMALL, mag(U())->weightedAverage(mesh_.V()).value());

    //Info<< "Uref - lRef: " << Uref << " - " << lRef_ << endl;

    scalar tau(lRef_/Uref_);

    tmp<volScalarField> B
    (
        CH3_ * Ostar() * max(Sstar(), Ostar())
        / max
        (
            (sqr(Sstar()) + sqr(Ostar())) / 2,
            dimensionedScalar
            (
                "Odt",
                dimensionSet(sqr(Sstar->dimensions())),
                sqr(OmegaLim_/tau)
            )
        )
    );


    tmp<volScalarField> g;

    if (approxTanh_)//faster
    {
        g = ptanh(B()*B()*B()*B());
    }
    else
    {
        B->dimensions().reset(dimless);
        g = tanh(B()*B()*B()*B());
    }
    g->dimensions().reset(dimless);

    B.clear();

    tmp<volScalarField> Kappa(sqrt(0.5*(sqr(Sstar()) + sqr(Ostar()))));
    Kappa.ref() = max
    (
        Kappa(),
        dimensionedScalar
        (
            "fscale",
            dimensionSet(Kappa->dimensions()),
            scalar(0.1)/tau
        )
    );
    Sstar.clear();
    Ostar.clear();

    scalar Cmu15 = Foam::pow(Cmu_, 1.5);

    tmp<volScalarField> lturb
    (
        sqrt(nuEff() / Kappa() / Cmu15)
    );
    Kappa.clear();

    lturb.ref().dimensions().reset(dimless);

    if (mesh_.time().outputTime() && debug)
    {
        //write all the constituent fields to file
        volScalarField averageg
        (
            IOobject
            (
                "averageg",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (fvc::average(g()))()
        );
        averageg.write();

        volScalarField CDESdelta
        (
            IOobject
            (
                "CDESdelta",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            CDES_ * delta()
        );
        CDESdelta.write();

        volScalarField averageLturb
        (
            IOobject
            (
                "averageLturb",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (fvc::average(lturb()))()
        );
        averageLturb.write();
    }

    tmp<volScalarField> ApowCH1
    (
        new volScalarField
        (
            IOobject
            (
                "ApowCH1",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("dl", dimensionSet(0,1,0,0,0,0,0), 0.0)
        )
    );

    ApowCH1.ref() = CDES_ * delta()
                  / max(lturb() * g(), 1e-15* dimensionedScalar("lr", lturb->dimensions(), lRef_))
                  - dimensionedScalar("o5", dimensionSet(0,1,0,0,0,0,0), 0.5);

    ApowCH1.ref().dimensions().reset(dimless);

    if (CH1_ == 3) //faster
    {
        ApowCH1.ref() = CH2_ * max(ApowCH1(), scalar(0.0));
        ApowCH1.ref() *= ApowCH1()*ApowCH1();
    }
    else
    {
        ApowCH1.ref() = Foam::pow(CH2_ * max(ApowCH1(), scalar(0.0)), CH1_);
    }


    if (approxTanh_)//faster
    {
        ApowCH1.ref() = ptanh(ApowCH1());
    }
    else
    {
        ApowCH1.ref() = tanh(ApowCH1());
    }

    tmp<surfaceScalarField> factor
    (
        new surfaceScalarField
        (
            IOobject
            (
                "blendingFactor",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("bf", dimless, 0.0),
            fixedValueFvsPatchScalarField::typeName
        )
    );

    ApowCH1.ref() = max(sigmaMax_* ApowCH1(), sigmaMin_);

    factor.ref() = fvc::interpolate(ApowCH1);

    return factor;
}

// ************************************************************************* //

} // End namespace blendingSources
} // End namespace Foam

// ************************************************************************* //
