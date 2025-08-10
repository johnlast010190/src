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

#include "submodels/MPPIC/ParticleStressModels/Lun/Lun.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ParticleStressModels
{
    defineTypeNameAndDebug(Lun, 0);

    addToRunTimeSelectionTable
    (
        ParticleStressModel,
        Lun,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ParticleStressModels::Lun::Lun
(
    const dictionary& dict
)
:
    ParticleStressModel(dict),
    e_(readScalar(dict.lookup("e"))),
    eps_(readScalar(dict.lookup("eps")))
{}


Foam::ParticleStressModels::Lun::Lun
(
    const Lun& ln
)
:
    ParticleStressModel(ln),
    e_(ln.e_),
    eps_(ln.eps_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ParticleStressModels::Lun::~Lun()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::ParticleStressModels::Lun::tau
(
    const Field<scalar>& alpha,
    const Field<scalar>& rho,
    const Field<scalar>& uSqr
) const
{
    tmp<Field<scalar>> g0
    (
        0.6
      / max
        (
            1.0 - pow(alpha/alphaPacked_,1.0/3.0),
            max(eps_*(1.0 - alpha), SMALL)
        )
    );

    tmp<Field<scalar>> gT(uSqr/3.0);

    return alpha*rho*(1.0 + alpha*(1.0 + e_)*g0)*gT;
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::ParticleStressModels::Lun::dTaudTheta
(
    const Field<scalar>& alpha,
    const Field<scalar>& rho,
    const Field<scalar>& uSqr
) const
{
    NotImplemented;
}


// ************************************************************************* //
