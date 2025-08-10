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
    (c) 2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/Kinematic/ParticleForces/Lift/TomiyamaLift/TomiyamaLiftForce.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::TomiyamaLiftForce<CloudType>::TomiyamaLiftForce::Cl
(
    const typename CloudType::parcelType& p,
    const vector& curlUc,
    const scalar Re,
    const scalar muc
) const
{
    const vector& g = this->owner().g().value();

    scalar Eo = p.Eo(g, p.d(), sigma_);
    scalar dH = p.d()*cbrt(1.0 + 0.163*pow(Eo, 0.757));
    scalar Eod = p.Eo(g, dH, sigma_);
    scalar f = 0.00105*pow3(Eod) - 0.0159*sqr(Eod) - 0.0204*Eod + 0.474;

    if (Eod <= 4)
    {
        return min(0.288*tanh(0.121*Re), f);
    }
    else if ((Eod > 4) && (Eod <= 10))
    {
        return f;
    }
    else
    {
        return -0.27;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TomiyamaLiftForce<CloudType>::TomiyamaLiftForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    LiftForce<CloudType>(owner, mesh, dict, forceType),
    sigma_(readScalar(this->coeffs().lookup("sigma")))
{}


template<class CloudType>
Foam::TomiyamaLiftForce<CloudType>::TomiyamaLiftForce
(
    const TomiyamaLiftForce& lf
)
:
    LiftForce<CloudType>(lf),
    sigma_(lf.sigma_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TomiyamaLiftForce<CloudType>::~TomiyamaLiftForce()
{}


// ************************************************************************* //
