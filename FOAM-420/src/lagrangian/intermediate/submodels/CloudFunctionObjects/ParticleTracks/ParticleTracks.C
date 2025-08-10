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

#include "submodels/CloudFunctionObjects/ParticleTracks/ParticleTracks.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "containers/Lists/ListListOps/ListListOps.H"
#include "db/IOobjects/IOPtrList/IOPtrList.H"

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleTracks<CloudType>::write()
{
    if (cloudPtr_.valid())
    {
        cloudPtr_->write();

        if (resetOnWrite_)
        {
            cloudPtr_->clear();
        }
    }
    else
    {
        if (debug)
        {
            InfoInFunction << "cloupPtr invalid" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleTracks<CloudType>::ParticleTracks
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    trackInterval_(readLabel(this->coeffDict().lookup("trackInterval"))),
    maxSamples_(readLabel(this->coeffDict().lookup("maxSamples"))),
    resetOnWrite_(this->coeffDict().lookup("resetOnWrite")),
    faceHitCounter_(),
    cloudPtr_(nullptr)
{}


template<class CloudType>
Foam::ParticleTracks<CloudType>::ParticleTracks
(
    const ParticleTracks<CloudType>& ppm
)
:
    CloudFunctionObject<CloudType>(ppm),
    trackInterval_(ppm.trackInterval_),
    maxSamples_(ppm.maxSamples_),
    resetOnWrite_(ppm.resetOnWrite_),
    faceHitCounter_(ppm.faceHitCounter_),
    cloudPtr_(ppm.cloudPtr_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleTracks<CloudType>::~ParticleTracks()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleTracks<CloudType>::preEvolve()
{
    if (!cloudPtr_.valid())
    {
        cloudPtr_.reset
        (
            this->owner().cloneBare(this->owner().name() + "Tracks").ptr()
        );
    }
}


template<class CloudType>
void Foam::ParticleTracks<CloudType>::postFace
(
    const parcelType& p,
    const label,
    bool&
)
{
    if
    (
        this->owner().solution().output()
     || this->owner().solution().transient()
    )
    {
        if (!cloudPtr_.valid())
        {
            FatalErrorInFunction
             << "Cloud storage not allocated" << abort(FatalError);
        }

        labelPairLookup::iterator iter =
            faceHitCounter_.find(labelPair(p.origProc(), p.origId()));

        label localI = -1;
        if (iter != faceHitCounter_.end())
        {
            iter()++;
            localI = iter();
        }
        else
        {
            localI = 1;
            faceHitCounter_.insert(labelPair(p.origProc(), p.origId()), localI);
        }

        label nSamples = floor(localI/trackInterval_);
        if ((localI % trackInterval_ == 0) && (nSamples < maxSamples_))
        {
            cloudPtr_->append
            (
                static_cast<parcelType*>(p.clone(this->owner().mesh()).ptr())
            );
        }
    }
}


// ************************************************************************* //
