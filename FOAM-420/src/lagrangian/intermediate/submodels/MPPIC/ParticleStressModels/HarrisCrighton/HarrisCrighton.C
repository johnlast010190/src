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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/MPPIC/ParticleStressModels/HarrisCrighton/HarrisCrighton.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ParticleStressModels
{
    defineTypeNameAndDebug(HarrisCrighton, 0);

    addToRunTimeSelectionTable
    (
        ParticleStressModel,
        HarrisCrighton,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ParticleStressModels::HarrisCrighton::HarrisCrighton
(
    const dictionary& dict
)
:
    ParticleStressModel(dict),
    pSolid_(readScalar(dict.lookup("pSolid"))),
    beta_(readScalar(dict.lookup("beta"))),
    eps_(readScalar(dict.lookup("eps")))
{}


Foam::ParticleStressModels::HarrisCrighton::HarrisCrighton
(
    const HarrisCrighton& hc
)
:
    ParticleStressModel(hc),
    pSolid_(hc.pSolid_),
    beta_(hc.beta_),
    eps_(hc.eps_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ParticleStressModels::HarrisCrighton::~HarrisCrighton()
{}


// * * * * * * * * * * * * * Privare Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::ParticleStressModels::HarrisCrighton::denominator
(
    const Field<scalar>& alpha
) const
{
    return
        max
        (
            alphaPacked_ - alpha,
            max(eps_*(1.0 - alpha), SMALL)
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::ParticleStressModels::HarrisCrighton::tau
(
    const Field<scalar>& alpha,
    const Field<scalar>& rho,
    const Field<scalar>& uSqr
) const
{
    return
    (
        pSolid_
      * pow(alpha, beta_)
      / denominator(alpha)
    );
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::ParticleStressModels::HarrisCrighton::dTaudTheta
(
    const Field<scalar>& alpha,
    const Field<scalar>& rho,
    const Field<scalar>& uSqr
) const
{
    const Field<scalar> d(denominator(alpha));

    return
    (
        pSolid_
      * pow(alpha, beta_)
      / d
      * (beta_/alpha + 1.0/d)
    );
}


// ************************************************************************* //
