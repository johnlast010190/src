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
    (c) 2016 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::CollidingParcel

Group
    grpLagrangianIntermediateParcels

Description
    Wrapper around kinematic parcel types to add collision modelling

SourceFiles
    CollidingParcelI.H
    CollidingParcel.C
    CollidingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef CollidingParcel_H
#define CollidingParcel_H

#include "particle/particle.H"
#include "parcels/Templates/CollidingParcel/CollisionRecordList/CollisionRecordList.H"
#include "fields/Fields/labelField/labelFieldIOField.H"
#include "fields/Fields/vectorField/vectorFieldIOField.H"
#include "primitives/demandDrivenEntry/demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef CollisionRecordList<vector, vector> collisionRecordList;
typedef vectorFieldCompactIOField pairDataFieldCompactIOField;
typedef vectorFieldCompactIOField wallDataFieldCompactIOField;

template<class ParcelType>
class CollidingParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const CollidingParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                       Class CollidingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class CollidingParcel
:
    public ParcelType
{
public:

    //- Size in bytes of the fields
    std::size_t sizeofFields() const;


    //- Class to hold thermo particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {

        // Private data

            //- Young's modulus [N/m2]
            demandDrivenEntry<scalar> youngsModulus_;

            //- Poisson's ratio
            demandDrivenEntry<scalar> poissonsRatio_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Member functions

            //- Return const access to Young's Modulus
            inline scalar youngsModulus() const;

            //- Return const access to Poisson's ratio
            inline scalar poissonsRatio() const;
    };


protected:

    // Protected data

        //- Force on particle due to collisions [N]
        vector f_;

        //- Angular momentum of Parcel in global reference frame [kg m2/s]
        vector angularMomentum_;

        //- Torque on particle due to collisions in global
        //  reference frame [Nm]
        vector torque_;

        //- Particle collision records
        collisionRecordList collisionRecords_;


public:

    // Static data members

        //- Runtime type information
        TypeName("CollidingParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " (fx fy fz)"
          + " (angularMomentumx angularMomentumy angularMomentumz)"
          + " (torquex torquey torquez)"
          + " collisionRecordsPairAccessed"
          + " collisionRecordsPairOrigProcOfOther"
          + " collisionRecordsPairOrigIdOfOther"
          + " (collisionRecordsPairData)"
          + " collisionRecordsWallAccessed"
          + " collisionRecordsWallPRel"
          + " (collisionRecordsWallData)"
        );

        //- String representation of property types
        AddToPropertyTypes
        (
            ParcelType,
            "{vector vector vector} *"
        );


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline CollidingParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from components
        inline CollidingParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const typename ParcelType::constantProperties& constProps
        );

        //- Construct from Istream
        CollidingParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        CollidingParcel(const CollidingParcel& p);

        //- Construct as a copy
        CollidingParcel(const CollidingParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new CollidingParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new CollidingParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<CollidingParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<CollidingParcel<ParcelType>>
                (
                    new CollidingParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to force
            inline const vector& f() const;

            //- Return const access to angular momentum
            inline const vector& angularMomentum() const;

            //- Return const access to torque
            inline const vector& torque() const;

            //- Return const access to the collision records
            inline const collisionRecordList& collisionRecords() const;

            //- Return access to force
            inline vector& f();

            //- Return access to angular momentum
            inline vector& angularMomentum();

            //- Return access to torque
            inline vector& torque();

            //- Return access to collision records
            inline collisionRecordList& collisionRecords();

            //- Particle angular velocity
            inline vector omega() const;


        // Tracking

            //- Move the parcel
            template<class TrackData>
            bool move(TrackData& td, const scalar trackTime);

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const transformer&);


       // I-O

            //- Read
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType>
            static void writeFields(const CloudType& c);

            //- Write particle fields as objects into the obr registry
            template<class CloudType>
            static void writeObjects(const CloudType& c, objectRegistry& obr);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const CollidingParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "parcels/Templates/CollidingParcel/CollidingParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "parcels/Templates/CollidingParcel/CollidingParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
