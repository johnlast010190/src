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
    (c) 2014-2015 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "incompressibleMultiphaseInteractingMixture.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleMultiphaseInteractingMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleMultiphaseInteractingMixture::
incompressibleMultiphaseInteractingMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    alphad_
    (
        IOobject
        (
            "alphad",
            U.time().timeName(),
            U.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("alphad", dimless, 0)
    ),
    phases_(lookup("phases"), phase::iNew(U, phi, this->lookupOrDefault("alphaMax", 1.0))),
    continuousPhase_(this->lookupOrDefault<word>("continuousPhase", "water")),
    gasPhase_(this->lookupOrDefault<word>("gasPhase", "none")),

    muModel_
    (
        mixtureViscosityModel::New
        (
            "mu",
            subDict("mixture"),
            U,
            phi
        )
    ),

    alphaMax_(this->lookupOrDefault("alphaMax", 1.0)),

    U_(U),
    phi_(phi),

    mu_
    (
        IOobject
        (
            "mu",
            U.time().timeName(),
            U.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),

    Hdm_
    (
        IOobject
        (
            "Hdm",
            U.time().timeName(),
            U.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedVector("Hdm", dimVelocity*dimDensity, vector::zero),
        calculatedFvPatchVectorField::typeName
    )
{
    if (U.mesh().foundObject<volScalarField>("alpha.air"))
    {
        gasPhase_ = "air";
    }
    Info<< "correcting mixture" << endl;
    correct();
    correctUdm();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::incompressibleMultiphaseInteractingMixture::correct()
{
    // compute disperse phases alpha
    alphad_.forceAssign(0.0);

    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        if (
            (iter().name() != continuousPhase_)
         && (iter().name() != gasPhase_)
        )
        {
            alphad_ += iter();
        }
    }

    // update material properties (viscosity)
    mu_ *= 0.0; // reset mu

    tmp<volScalarField> tmuc
    (
        new volScalarField
        (
            IOobject
            (
                "muc",
                alphad_.time().timeName(),
                alphad_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            alphad_.mesh(),
            dimensionedScalar("muc", dimensionSet(1, -1, -1, 0, 0), 0),
            calculatedFvPatchScalarField::typeName
        )
    );
    volScalarField& muc = tmuc.ref();
    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        if (iter().name() == continuousPhase_)
        {
            muc += iter().rho()*iter().nu();
        }
    }

    mu_ = muModel_->mu(tmuc);

    // linear blending with air viscosity (if air phase present)
    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        if (iter().name() == gasPhase_)
        {
            mu_ *= (1.-iter());
            mu_ += iter()*iter().rho()*iter().nu();
        }
    }
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleMultiphaseInteractingMixture::rho() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> trho = iter()*iter().rho();
    volScalarField& rho = trho.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        rho += iter()*iter().rho();
    }

    return trho;
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleMultiphaseInteractingMixture::rho(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> trho = iter().boundaryField()[patchi]*iter().rho().value();
    scalarField& rho = trho.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        rho += iter().boundaryField()[patchi]*iter().rho().value();
    }

    return trho;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::incompressibleMultiphaseInteractingMixture::tauDrift() const
{
    tmp<volVectorField> Hm = Hmass();
    tmp<volVectorField> H = Hvol();

    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volSymmTensorField> ttauDrift = iter()*iter().rho()*sqr(iter().Urel());
    volSymmTensorField& tauDrift = ttauDrift.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        tauDrift += iter()*iter().rho()*sqr(iter().Urel());
    }

    tauDrift += rho()*sqr(H())
      + rho()*(sqr(H()+U()) - sqr(U()) - sqr(H()))
      - rho()*(sqr(U()+(Hm()/rho())) - sqr(Hm()/rho()) - sqr(U()))
      - rho()*(sqr((Hm()/rho())+H()) - sqr(Hm()/rho()) - sqr(H()));

    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            "tauDm",
            ttauDrift
        )
    );
}


Foam::tmp<Foam::volVectorField>
Foam::incompressibleMultiphaseInteractingMixture::Hmass() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volVectorField> tHmass = iter()*iter().rho()*iter().Urel();
    volVectorField& Hmass = tHmass.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Hmass += iter()*iter().rho()*iter().Urel();
    }

    return tHmass;
}


Foam::tmp<Foam::volVectorField>
Foam::incompressibleMultiphaseInteractingMixture::Hvol() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volVectorField> tHvol = iter()*iter().Urel();
    volVectorField& Hvol = tHvol.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Hvol += iter()*iter().Urel();
    }

    return tHvol;
}


bool Foam::incompressibleMultiphaseInteractingMixture::read()
{
    if (regIOobject::read())
    {
        alphaMax_ =
            this->lookupOrDefault
            (
                "alphaMax",
                1.0
            );

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
