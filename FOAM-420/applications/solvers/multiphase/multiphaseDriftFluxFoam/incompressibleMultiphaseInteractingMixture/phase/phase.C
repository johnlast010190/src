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
    (c) 2010-2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "phase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phase::phase
(
    const word& phaseName,
    const dictionary& phaseDict,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const scalar& alphaMax
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    name_(phaseName),
    phaseDict_(phaseDict),
    alphaMax_(alphaMax),
    useTotalSolids_(phaseDict.lookupOrDefault<Switch>("useTotalSolids", false)),
    rho_(nullptr),
    nu_(nullptr),
    UdmModelPtr_()
{
    read(phaseDict);

    if (useTotalSolids_)
    {
        UdmModelPtr_ =
            relativeVelocityModel::New
            (
                name_,
                phaseDict_,
                this->mesh().lookupObject<volScalarField>("alphad"),
                alphaMax_
            );
    }
    else
    {
        UdmModelPtr_ =
            relativeVelocityModel::New
            (
                name_,
                phaseDict_,
                *this,
                alphaMax_
            );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phase> Foam::phase::clone() const
{
    NotImplemented;
}


void Foam::phase::correctUdm()
{
    UdmModelPtr_->correct();
}


bool Foam::phase::read(const dictionary& phaseDict)
{
    phaseDict_ = phaseDict;

    rho_.reset
    (
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rho",
                phaseDict_,
                dimDensity,
                1000.0
            )
        )
    );

    nu_.reset
    (
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "nu",
                phaseDict_,
                dimArea/dimTime,
                1e-06
            )
        )
    );

    return true;
}


// ************************************************************************* //
