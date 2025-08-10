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

#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "viscosityModels/viscosityModel/viscosityModel.H"
#include "fields/volFields/volFields.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singlePhaseTransportModel, 0);
}


const Foam::physicalProperties&
Foam::singlePhaseTransportModel::properties() const
{
    return physicalPropertiesPtr_();
}

//- eturn const access to viscosityModel
const Foam::viscosityModel& Foam::singlePhaseTransportModel::viscosity() const
{
    return viscosityModelPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singlePhaseTransportModel::singlePhaseTransportModel
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
    viscosityModelPtr_(viscosityModel::New("nu", *this, U, phi)),
    physicalPropertiesPtr_(new physicalProperties(*this))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singlePhaseTransportModel::~singlePhaseTransportModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::singlePhaseTransportModel::nu() const
{
    return viscosity().nu();
}

Foam::tmp<Foam::scalarField>
Foam::singlePhaseTransportModel::nu(const label patchi) const
{
    return viscosity().nu(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::singlePhaseTransportModel::rho() const
{
    dimensionedScalar rho = properties().rho();

    tmp<volScalarField> rhot
    (
        new volScalarField
        (
            IOobject
            (
                rho.name(),
                nu()().time().timeName(),
                nu()().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            nu()().mesh(),
            rho,
            extrapolatedCalculatedFvPatchField<scalar>::typeName
        )
    );

    return rhot;
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseTransportModel::Cp() const
{
    dimensionedScalar Cp = properties().Cp();

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                Cp.name(),
                nu()().time().timeName(),
                nu()().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            nu()().mesh(),
            Cp,
            extrapolatedCalculatedFvPatchField<scalar>::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseTransportModel::lambda() const
{
    dimensionedScalar lambda = properties().lambda();

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                lambda.name(),
                nu()().time().timeName(),
                nu()().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            nu()().mesh(),
            lambda,
            extrapolatedCalculatedFvPatchField<scalar>::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseTransportModel::Prt() const
{
    dimensionedScalar Prt = properties().Prt();

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                Prt.name(),
                nu()().time().timeName(),
                nu()().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            nu()().mesh(),
            Prt,
            extrapolatedCalculatedFvPatchField<scalar>::typeName
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::singlePhaseTransportModel::Pr() const
{

    if (properties().PrValid())
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Pr",
                    nu()().time().timeName(),
                    nu()().db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                nu()().mesh(),
                properties().Pr(),
                extrapolatedCalculatedFvPatchField<scalar>::typeName
            )
        );
    }
    else if
    (
        properties().CpValid()
     && properties().rhoValid()
     && properties().lambdaValid()
    )
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Pr",
                    nu()().time().timeName(),
                    nu()().db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                nu()()*properties().rho()
                *properties().Cp()/properties().lambda(),
                extrapolatedCalculatedFvPatchField<scalar>::typeName
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Constants necessary for Prandtl number calculation"
            << "not defined in transportProperties." << nl
            << "Requires either k/rho/Cp/nu or Pr."
            << exit(FatalError);
    }

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Pr",
                nu()().time().timeName(),
                nu()().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
                nu()().mesh(),
                properties().Pr(),
            extrapolatedCalculatedFvPatchField<scalar>::typeName
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::singlePhaseTransportModel::alphaLam() const
{

    if
    (
        properties().CpValid()
     && properties().rhoValid()
     && properties().lambdaValid()
    )
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "alpha",
                    nu()().time().timeName(),
                    nu()().db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                nu()().mesh(),
                properties().lambda()/properties().rho()/properties().Cp(),
                extrapolatedCalculatedFvPatchField<scalar>::typeName
            )
        );
    }
    else if (properties().PrValid())
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "alpha",
                    nu()().time().timeName(),
                    nu()().db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                nu()()/properties().Pr(),
                extrapolatedCalculatedFvPatchField<scalar>::typeName
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Constants necessary for thermal diffusivity calculation"
            << "not defined in transportProperties." << nl
            << "Requires either k/rho/Cp or nu/Pr."
            << exit(FatalError);
    }

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "alpha",
                nu()().time().timeName(),
                nu()().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            nu()().mesh(),
            properties().lambda()/properties().rho()/properties().Cp(),
            extrapolatedCalculatedFvPatchField<scalar>::typeName
        )
    );

}


void Foam::singlePhaseTransportModel::correct()
{
    viscosityModelPtr_->correct();
    physicalPropertiesPtr_->correct();
}


bool Foam::singlePhaseTransportModel::read()
{
    if (regIOobject::read())
    {
        return
        (
            viscosityModelPtr_->read(*this)
            && physicalPropertiesPtr_->read(*this)
        );
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
