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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::MPPICParcel<ParcelType>::MPPICParcel
(
    const polyMesh& owner,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(owner, position, celli, tetFacei, tetPti),
    UCorrect_(Zero)
{}


template<class ParcelType>
inline Foam::MPPICParcel<ParcelType>::MPPICParcel
(
    const polyMesh& owner,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& UCorrect0,
    const typename ParcelType::constantProperties& constProps
)
:
    ParcelType
    (
        owner,
        position,
        celli,
        tetFacei,
        tetPti,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        constProps
    ),
    UCorrect_(UCorrect0)
{}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::vector& Foam::MPPICParcel<ParcelType>::UCorrect() const
{
    return UCorrect_;
}


template<class ParcelType>
inline Foam::vector& Foam::MPPICParcel<ParcelType>::UCorrect()
{
    return UCorrect_;
}


// ************************************************************************* //
