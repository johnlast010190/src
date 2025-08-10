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
    (c) 2012-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/Kinematic/ParticleForces/Lift/SaffmanMeiLift/SaffmanMeiLiftForce.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::SaffmanMeiLiftForce<CloudType>::SaffmanMeiLiftForce::Cl
(
    const typename CloudType::parcelType& p,
    const vector& curlUc,
    const scalar Re,
    const scalar muc
) const
{
    scalar Rew = p.rhoc()*mag(curlUc)*sqr(p.d())/(muc + ROOTVSMALL);
    scalar beta = 0.5*(Rew/(Re + ROOTVSMALL));
    scalar alpha = 0.3314*sqrt(beta);
    scalar f = (1.0 - alpha)*exp(-0.1*Re) + alpha;

    scalar Cld = 0.0;
    if (Re < 40)
    {
        Cld = 6.46*f;
    }
    else
    {
        Cld = 6.46*0.0524*sqrt(beta*Re);
    }

    return 3.0/(mathematical::twoPi*sqrt(Rew + ROOTVSMALL))*Cld;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SaffmanMeiLiftForce<CloudType>::SaffmanMeiLiftForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    LiftForce<CloudType>(owner, mesh, dict, forceType)
{}


template<class CloudType>
Foam::SaffmanMeiLiftForce<CloudType>::SaffmanMeiLiftForce
(
    const SaffmanMeiLiftForce& lf
)
:
    LiftForce<CloudType>(lf)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SaffmanMeiLiftForce<CloudType>::~SaffmanMeiLiftForce()
{}


// ************************************************************************* //
