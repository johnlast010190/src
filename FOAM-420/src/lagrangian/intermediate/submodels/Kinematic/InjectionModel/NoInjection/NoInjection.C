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

\*---------------------------------------------------------------------------*/

#include "submodels/Kinematic/InjectionModel/NoInjection/NoInjection.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoInjection<CloudType>::NoInjection
(
    const dictionary&,
    CloudType& owner,
    const word&
)
:
    InjectionModel<CloudType>(owner)
{}


template<class CloudType>
Foam::NoInjection<CloudType>::NoInjection(const NoInjection<CloudType>& im)
:
    InjectionModel<CloudType>(im.owner_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoInjection<CloudType>::~NoInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NoInjection<CloudType>::active() const
{
    return false;
}


template<class CloudType>
Foam::scalar Foam::NoInjection<CloudType>::timeEnd() const
{
    return 0.0;
}


template<class CloudType>
Foam::label Foam::NoInjection<CloudType>::parcelsToInject
(
    const scalar,
    const scalar
)
{
    return 0;
}


template<class CloudType>
Foam::scalar Foam::NoInjection<CloudType>::volumeToInject
(
    const scalar,
    const scalar
)
{
    return 0.0;
}


template<class CloudType>
void Foam::NoInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector&,
    label&,
    label&,
    label&
)
{}


template<class CloudType>
void Foam::NoInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() = Zero;

    // set particle diameter
    parcel.d() = 0.0;
}


template<class CloudType>
bool Foam::NoInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::NoInjection<CloudType>::validInjection(const label)
{
    return false;
}


// ************************************************************************* //
