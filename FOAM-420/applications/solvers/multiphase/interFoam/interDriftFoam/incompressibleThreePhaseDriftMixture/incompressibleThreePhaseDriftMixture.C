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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2010-2012 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "incompressibleThreePhaseDriftMixture.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvc.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::incompressibleThreePhaseDriftMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();
    nuModel3_->correct();

    volScalarField limitedAlpha1(max(alpha1_, scalar(0)));
    volScalarField limitedAlpha2(max(alpha2_, scalar(0)));
    volScalarField limitedAlpha3(max(alpha3_, scalar(0)));

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(limitedAlpha1*physProp1_->rho() + limitedAlpha2*physProp2_->rho()
                + limitedAlpha3*physProp3_->rho());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleThreePhaseDriftMixture::incompressibleThreePhaseDriftMixture
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

    phase1Name_(wordList(lookup("phases"))[0]),
    phase2Name_(wordList(lookup("phases"))[1]),
    phase3Name_(wordList(lookup("phases"))[2]),

    alpha1_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase1Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha3_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase3Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    U_(U),
    phi_(phi),

    nu_
    (
        IOobject
        (
            "nu",
            U.time().timeName(),
            U.db()
        ),
        U.mesh(),
        dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
    nuModel3_
    (
        viscosityModel::New
        (
            "nu3",
            subDict(phase3Name_),
            U,
            phi
        )
    ),

    physProp1_(new physicalProperties(subDict(phase1Name_))),
    physProp2_(new physicalProperties(subDict(phase2Name_))),
    physProp3_(new physicalProperties(subDict(phase3Name_))),

    UdmModelPtr_
    (
        relativeVelocityModel::New
        (
            "driftVel2",
            subDict(phase2Name_),
            alpha2_
        )
    )
{
    alpha3_.forceAssign(1.0 - alpha1_ - alpha2_);
    calcNu();
    correctUdm();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseDriftMixture::mu() const
{
    volScalarField limitedAlpha1(max(alpha1_, scalar(0)));
    volScalarField limitedAlpha2(max(alpha2_, scalar(0)));
    volScalarField limitedAlpha3(max(alpha3_, scalar(0)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            limitedAlpha1*physProp1_->rho()*nuModel1_->nu()
          + limitedAlpha2*physProp2_->rho()*nuModel2_->nu()
          + limitedAlpha3*physProp3_->rho()*nuModel3_->nu()
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseDriftMixture::muf() const
{
    volScalarField limitedAlpha1(max(alpha1_, scalar(0)));
    volScalarField limitedAlpha2(max(alpha2_, scalar(0)));
    volScalarField limitedAlpha3(max(alpha3_, scalar(0)));

    surfaceScalarField alpha1f(fvc::interpolate(limitedAlpha1));
    surfaceScalarField alpha2f(fvc::interpolate(limitedAlpha2));
    surfaceScalarField alpha3f(fvc::interpolate(limitedAlpha3));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "mu",
            alpha1f*physProp1_->rho()*fvc::interpolate(nuModel1_->nu())
          + alpha2f*physProp2_->rho()*fvc::interpolate(nuModel2_->nu())
          + alpha3f*physProp3_->rho()*fvc::interpolate(nuModel3_->nu())
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseDriftMixture::nuf() const
{
    volScalarField limitedAlpha1(max(alpha1_, scalar(0)));
    volScalarField limitedAlpha2(max(alpha2_, scalar(0)));
    volScalarField limitedAlpha3(max(alpha3_, scalar(0)));

    surfaceScalarField alpha1f(fvc::interpolate(limitedAlpha1));
    surfaceScalarField alpha2f(fvc::interpolate(limitedAlpha2));
    surfaceScalarField alpha3f(fvc::interpolate(limitedAlpha3));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nu",
            (
                alpha1f*physProp1_->rho()*fvc::interpolate(nuModel1_->nu())
              + alpha2f*physProp2_->rho()*fvc::interpolate(nuModel2_->nu())
              + alpha3f*physProp3_->rho()*fvc::interpolate(nuModel3_->nu())
            )/(alpha1f*physProp1_->rho() + alpha2f*physProp2_->rho()
               + alpha3f*physProp3_->rho())
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseDriftMixture::rho() const
{
    volScalarField limitedAlpha1(max(alpha1_, scalar(0)));
    volScalarField limitedAlpha2(max(alpha2_, scalar(0)));
    volScalarField limitedAlpha3(max(alpha3_, scalar(0)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "rho",
            limitedAlpha1*physProp1_->rho()
          + limitedAlpha2*physProp2_->rho()
          + limitedAlpha3*physProp3_->rho()
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseDriftMixture::Cp() const
{
    NotImplemented;

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            alpha1_*physProp1_->rho()*physProp1_->Cp()
          + alpha2_*physProp2_->rho()*physProp2_->Cp()
          + alpha3_*physProp3_->rho()*physProp3_->Cp()
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseDriftMixture::lambda() const
{
    NotImplemented;

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "lambda",
            alpha1_*physProp1_->rho()*physProp1_->lambda()
          + alpha2_*physProp2_->rho()*physProp2_->lambda()
          + alpha3_*physProp3_->rho()*physProp3_->lambda()
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseDriftMixture::Prt() const
{
    NotImplemented;

    volScalarField limitedAlpha1( min(max(alpha1_, scalar(0)), scalar(1)) );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "Prt",
            alpha1_*physProp1_->rho()*physProp1_->Prt()
          + alpha2_*physProp2_->rho()*physProp2_->Prt()
          + alpha3_*physProp3_->rho()*physProp3_->Prt()
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseDriftMixture::Pr() const
{
    NotImplemented;

    volScalarField limitedAlpha1( min(max(alpha1_, scalar(0)), scalar(1)) );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "Pr",
            nu()/alphaLam()
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseDriftMixture::alphaLam() const
{
    NotImplemented;

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "alpha",
            lambda()/rho()/Cp()
        )
    );
}



bool Foam::incompressibleThreePhaseDriftMixture::read()
{
    if (transportModel::read())
    {
        if
        (
            nuModel1_().read(subDict(phase1Name_))
         && nuModel2_().read(subDict(phase2Name_))
         && nuModel3_().read(subDict(phase3Name_))
         && physProp1_().read(subDict(phase1Name_))
         && physProp2_().read(subDict(phase2Name_))
         && physProp3_().read(subDict(phase3Name_))
        )
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
