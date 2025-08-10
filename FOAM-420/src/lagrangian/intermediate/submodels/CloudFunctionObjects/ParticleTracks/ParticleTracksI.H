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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/CloudFunctionObjects/ParticleTracks/ParticleTracks.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
inline Foam::label Foam::ParticleTracks<CloudType>::trackInterval() const
{
    return trackInterval_;
}


template<class CloudType>
inline Foam::label Foam::ParticleTracks<CloudType>::maxSamples() const
{
    return maxSamples_;
}


template<class CloudType>
inline const Foam::Switch& Foam::ParticleTracks<CloudType>::resetOnWrite() const
{
    return resetOnWrite_;
}


template<class CloudType>
inline const Foam::labelPairLookup&
Foam::ParticleTracks<CloudType>::faceHitCounter() const
{
    return faceHitCounter_;
}


template<class CloudType>
inline const Foam::Cloud<typename CloudType::parcelType>&
Foam::ParticleTracks<CloudType>::cloud() const
{
    return cloudPtr_();
}


// ************************************************************************* //
