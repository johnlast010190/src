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

#include "incompressibleThreePhaseMixture.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvc.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::incompressibleThreePhaseMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();
    nuModel3_->correct();

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(alpha1lim()*physProp1_->rho() + alpha2lim()*physProp2_->rho()
                + alpha3lim()*physProp3_->rho());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleThreePhaseMixture::incompressibleThreePhaseMixture
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
    physProp3_(new physicalProperties(subDict(phase3Name_)))
{
    alpha3_.forceAssign(1.0 - alpha1_ - alpha2_);
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::mu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            alpha1lim()*physProp1_->rho()*nuModel1_->nu()
          + alpha2lim()*physProp2_->rho()*nuModel2_->nu()
          + alpha3lim()*physProp3_->rho()*nuModel3_->nu()
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::muf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));

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
Foam::incompressibleThreePhaseMixture::nuf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));

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
Foam::incompressibleThreePhaseMixture::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "rho",
            alpha1lim()*physProp1_->rho()
          + alpha2lim()*physProp2_->rho()
          + alpha3lim()*physProp3_->rho()
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::Cp() const
{
    NotImplemented;
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::lambda() const
{
    NotImplemented;
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::Prt() const
{
    NotImplemented;
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::Pr() const
{
    NotImplemented;
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::alphaLam() const
{
    NotImplemented;
}



bool Foam::incompressibleThreePhaseMixture::read()
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
