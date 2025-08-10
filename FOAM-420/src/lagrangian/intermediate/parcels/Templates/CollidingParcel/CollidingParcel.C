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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "parcels/Templates/CollidingParcel/CollidingParcel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::CollidingParcel<ParcelType>::CollidingParcel
(
    const CollidingParcel<ParcelType>& p
)
:
    ParcelType(p),
    f_(p.f_),
    angularMomentum_(p.angularMomentum_),
    torque_(p.torque_),
    collisionRecords_(p.collisionRecords_)
{}


template<class ParcelType>
Foam::CollidingParcel<ParcelType>::CollidingParcel
(
    const CollidingParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    f_(p.f_),
    angularMomentum_(p.angularMomentum_),
    torque_(p.torque_),
    collisionRecords_(p.collisionRecords_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::CollidingParcel<ParcelType>::move
(
    TrackData& td,
    const scalar trackTime
)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    td.keepParticle = true;
    td.switchProcessor = false;

    switch (td.part())
    {
        case TrackData::tpVelocityHalfStep:
        {
            // First and last leapfrog velocity adjust part, required
            // before and after tracking and force calculation

            p.U() += 0.5*trackTime*p.f()/p.mass();

            p.angularMomentum() += 0.5*trackTime*p.torque();

            td.keepParticle = true;

            p.face() = -1;

            break;
        }

        case TrackData::tpLinearTrack:
        {
            ParcelType::move(td, trackTime);

            break;
        }

        case TrackData::tpRotationalTrack:
        {
            NotImplemented;

            break;
        }

        default:
        {
            FatalErrorInFunction
                << td.part() << " is an invalid part of the tracking method."
                << abort(FatalError);
        }
    }

    return td.keepParticle;
}


template<class ParcelType>
void Foam::CollidingParcel<ParcelType>::transformProperties
(
    const transformer& transform
)
{
    ParcelType::transformProperties(transform);

    f_ = transform.transform(f_);
    angularMomentum_ = transform.transform(angularMomentum_);
    torque_ = transform.transform(torque_);
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "parcels/Templates/CollidingParcel/CollidingParcelIO.C"

// ************************************************************************* //
