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
    (c) 2010-2012 Esi Ltd.
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/Kinematic/ParticleForces/ElectricField/ElectricFieldForce.H"
#include "include/demandDrivenData.H"
#include "global/constants/electromagnetic/electromagneticConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ElectricFieldForce<CloudType>::ElectricFieldForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    EFieldName_
    (
        this->coeffs().template lookupOrDefault<word>("EField", "E")
    ),
    EFieldInterpPtr_(nullptr),
    particleCharge_
    (
        readScalar(this->coeffs().lookup("particleCharge"))
    )
{}


template<class CloudType>
Foam::ElectricFieldForce<CloudType>::ElectricFieldForce
(
    const ElectricFieldForce& pf
)
:
    ParticleForce<CloudType>(pf),
    EFieldName_(pf.EFieldName_),
    EFieldInterpPtr_(pf.EFieldInterpPtr_),
    particleCharge_(pf.particleCharge_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ElectricFieldForce<CloudType>::~ElectricFieldForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ElectricFieldForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const volVectorField& EField =
            this->mesh().template lookupObject<volVectorField>(EFieldName_);

        EFieldInterpPtr_ = interpolation<vector>::New
        (
            this->owner().solution().interpolationSchemes(),
            EField
        ).ptr();
    }
    else
    {
        deleteDemandDrivenData(EFieldInterpPtr_);
    }
}


template<class CloudType>
Foam::forceSuSp Foam::ElectricFieldForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(vector::zero, 0.0);

    const interpolation<vector>& EFieldInterp = *EFieldInterpPtr_;

    value.Su()= particleCharge_ * p.nParticle()
                *EFieldInterp.interpolate(p.position(), p.currentTetIndices());

    return value;
}


// ************************************************************************* //
