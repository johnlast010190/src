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

Class
    Foam::ParticleTracks

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Records particle state (all variables) on each call to postFace

SourceFiles
    ParticleTracks.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleTracks_H
#define ParticleTracks_H

#include "submodels/CloudFunctionObjects/CloudFunctionObject/CloudFunctionObject.H"
#include "primitives/Pair/labelPairHashes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleTracks Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleTracks
:
    public CloudFunctionObject<CloudType>
{
    // Private data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Number of face-hit intervals between storing parcel data
        label trackInterval_;

        //- Maximum number of particles to store per track
        label maxSamples_;

        //- Flag to indicate whether data should be reset/cleared on writing
        Switch resetOnWrite_;

        //- Table of number of times a particle has hit a face
        labelPairLookup faceHitCounter_;

        //- Pointer to the cloud storage
        autoPtr<Cloud<parcelType>> cloudPtr_;


protected:

    // Protected member functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("particleTracks");


    // Constructors

        //- Construct from dictionary
        ParticleTracks
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParticleTracks(const ParticleTracks<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new ParticleTracks<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleTracks();


    // Member Functions

        // Access

            //- Return const access to the track interval
            inline label trackInterval() const;

            //- Return const access to the max samples
            inline label maxSamples() const;

            //- Return const access to the reset on write flag
            inline const Switch& resetOnWrite() const;

            //- Return the table of number of times a particle has hit a face
            inline const labelPairLookup& faceHitCounter() const;

            //- Return const access to the cloud
            inline const Cloud<parcelType>& cloud() const;


        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Post-face hook
            virtual void postFace
            (
                const parcelType& p,
                const label facei,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/CloudFunctionObjects/ParticleTracks/ParticleTracksI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/CloudFunctionObjects/ParticleTracks/ParticleTracks.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
