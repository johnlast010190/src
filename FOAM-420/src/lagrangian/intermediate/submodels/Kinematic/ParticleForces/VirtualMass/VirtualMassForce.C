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

#include "submodels/Kinematic/ParticleForces/VirtualMass/VirtualMassForce.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VirtualMassForce<CloudType>::VirtualMassForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    PressureGradientForce<CloudType>(owner, mesh, dict, forceType),
    Cvm_(readScalar(this->coeffs().lookup("Cvm")))
{}


template<class CloudType>
Foam::VirtualMassForce<CloudType>::VirtualMassForce
(
    const VirtualMassForce& vmf
)
:
    PressureGradientForce<CloudType>(vmf),
    Cvm_(vmf.Cvm_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VirtualMassForce<CloudType>::~VirtualMassForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::VirtualMassForce<CloudType>::cacheFields(const bool store)
{
    PressureGradientForce<CloudType>::cacheFields(store);
}


template<class CloudType>
Foam::forceSuSp Foam::VirtualMassForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value =
        PressureGradientForce<CloudType>::calcCoupled(p, dt, mass, Re, muc);

    value.Su() *= Cvm_;

    return value;
}


template<class CloudType>
Foam::scalar Foam::VirtualMassForce<CloudType>::massAdd
(
    const typename CloudType::parcelType& p,
    const scalar mass
) const
{
    return mass*p.rhoc()/p.rho()*Cvm_;
}


// ************************************************************************* //
