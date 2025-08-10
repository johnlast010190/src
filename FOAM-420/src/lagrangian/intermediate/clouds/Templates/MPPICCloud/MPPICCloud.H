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

Class
    Foam::MPPICCloud

Group
    grpLagrangianIntermediateClouds

Description
    Adds MPPIC modelling to kinematic clouds

SourceFiles
    MPPICCloudI.H
    MPPICCloud.C

\*---------------------------------------------------------------------------*/

#ifndef MPPICCloud_H
#define MPPICCloud_H

#include "particle/particle.H"
#include "Cloud/Cloud.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class PackingModel;

template<class CloudType>
class DampingModel;

template<class CloudType>
class IsotropyModel;

/*---------------------------------------------------------------------------*\
                       Class MPPICCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MPPICCloud
:
    public CloudType
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::parcelType parcelType;

        //- Convenience typedef for this cloud type
        typedef MPPICCloud<CloudType> MPPICCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<MPPICCloud<CloudType>> cloudCopyPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        MPPICCloud(const MPPICCloud&);

        //- Disallow default bitwise assignment
        void operator=(const MPPICCloud&);


protected:

    // Protected data

        // References to the cloud sub-models

            //- Packing model
            autoPtr<PackingModel<MPPICCloud<CloudType>>> packingModel_;

            //- Damping model
            autoPtr<DampingModel<MPPICCloud<CloudType>>>
                dampingModel_;

            //- Exchange model
            autoPtr<IsotropyModel<MPPICCloud<CloudType>>>
                isotropyModel_;


        // Initialisation

            //- Set cloud sub-models
            void setModels();


public:

    // Constructors

        //- Construct given carrier gas fields
        MPPICCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        MPPICCloud
        (
            MPPICCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        MPPICCloud
        (
            const fvMesh& mesh,
            const word& name,
            const MPPICCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new MPPICCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new MPPICCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~MPPICCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const MPPICCloud& cloudCopy() const;

            //- Return const access to the packing model
            inline const PackingModel<MPPICCloud<CloudType>>&
                packingModel() const;

            //- Return a reference to the packing model
            inline PackingModel<MPPICCloud<CloudType>>& packingModel();

            //- Return condt access to the damping model
            inline const DampingModel<MPPICCloud<CloudType>>&
                dampingModel() const;

            //- Return a reference to the damping model
            inline DampingModel<MPPICCloud<CloudType>>& dampingModel();

            //- Return condt access to the isotropy model
            inline const IsotropyModel<MPPICCloud<CloudType>>&
                isotropyModel() const;

            //- Return a reference to the isotropy model
            inline IsotropyModel<MPPICCloud<CloudType>>& isotropyModel();


        // Cloud evolution functions

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Evolve the cloud
            void evolve();

            //- Particle motion
            template<class TrackData>
            void motion(TrackData& td);


        //- I-O

            //- Print cloud information
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "clouds/Templates/MPPICCloud/MPPICCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "clouds/Templates/MPPICCloud/MPPICCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
