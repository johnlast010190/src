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
    Foam::WallCollisionRecord

Description
    Record of a collision between the particle holding the record and
    a wall face at the position relative to the centre of the particle.

SourceFiles
    WallCollisionRecordI.H
    WallCollisionRecord.C
    WallCollisionRecordIO.C

\*---------------------------------------------------------------------------*/

#ifndef WallCollisionRecord_H
#define WallCollisionRecord_H

#include "primitives/Vector/vector/vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type>
class WallCollisionRecord;

template<class Type>
inline bool operator==
(
    const WallCollisionRecord<Type>&,
    const WallCollisionRecord<Type>&
);

template<class Type>
inline bool operator!=
(
    const WallCollisionRecord<Type>&,
    const WallCollisionRecord<Type>&
);

template<class Type>
Istream& operator>>(Istream&, WallCollisionRecord<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const WallCollisionRecord<Type>&);


/*---------------------------------------------------------------------------*\
                         Class WallCollisionRecord Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class WallCollisionRecord
{
    // Private data

        // //- Recording whether or not this record has been accessed
        bool accessed_;

        //- The position of wall impact relative to the particle centre
        vector pRel_;

        //- Collision data, stored as if the storing particle was the
        //  first particle (particle A) in the collision.
        Type data_;


public:

    // Static data members

        //- Tolerance for detecting seriously erroneous wall matches
        static const scalar errorCosAngle;


    // Constructors

        //- Construct null
        WallCollisionRecord();

        //- Construct from components
        WallCollisionRecord
        (
            bool accessed,
            const vector& pRel,
            const Type& data = Type(Zero)
        );

        //- Construct from Istream
        WallCollisionRecord(Istream&);

        //- Construct as copy
        WallCollisionRecord(const WallCollisionRecord&);


    //- Destructor
    ~WallCollisionRecord();


    // Member Functions


        // Access

            //- Return the pRel data
            inline const vector& pRel() const;

            //- Return access to the collision data
            inline const Type& collisionData() const;

            //- Return access to the collision data
            inline Type& collisionData();


        // Check

            inline bool match(const vector& pRel, scalar radius);

            //- Return the accessed status of the record
            inline bool accessed() const;


        // Edit

            //- Set the accessed property of the record to accessed
            inline void setAccessed();

            //- Set the accessed property of the record to unaccessed
            inline void setUnaccessed();


    // Member Operators

        void operator=(const WallCollisionRecord&);


    // Friend Operators

        friend bool operator== <Type>
        (
            const WallCollisionRecord<Type>&,
            const WallCollisionRecord<Type>&
        );

        friend bool operator!= <Type>
        (
            const WallCollisionRecord<Type>&,
            const WallCollisionRecord<Type>&
        );


    // IOstream Operators

        friend Istream& operator>> <Type>
        (
            Istream&,
            WallCollisionRecord<Type>&
        );

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const WallCollisionRecord<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "parcels/Templates/CollidingParcel/CollisionRecordList/WallCollisionRecord/WallCollisionRecordI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "parcels/Templates/CollidingParcel/CollisionRecordList/WallCollisionRecord/WallCollisionRecord.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
