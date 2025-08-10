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
    Foam::TrajectoryCollision

Description
    Trajectory collision model by N. Nordin, based on O'Rourke's collision
    model

\*---------------------------------------------------------------------------*/

#ifndef TrajectoryCollision_H
#define TrajectoryCollision_H

#include "submodels/StochasticCollision/ORourkeCollision/ORourkeCollision.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class TrajectoryCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TrajectoryCollision
:
    public ORourkeCollision<CloudType>
{
protected:

    // Protected Data

        //- Convenience typedef to the cloud's parcel type
        typedef typename CloudType::parcelType parcelType;

        //- Space coefficient
        scalar cSpace_;

        //- Time coefficient
        scalar cTime_;


    // Protected Member Functions

        //- Main collision routine
        virtual void collide(const scalar dt);

        //- Collide parcels and return true if mass has changed
        virtual bool collideParcels
        (
            const scalar dt,
            parcelType& p1,
            parcelType& p2,
            scalar& m1,
            scalar& m2
        );


public:

    //- Runtime type information
    TypeName("trajectory");


    // Constructors

        //- Construct from dictionary
        TrajectoryCollision(const dictionary& dict, CloudType& cloud);

        //- Construct copy
        TrajectoryCollision(const TrajectoryCollision<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<StochasticCollisionModel<CloudType>> clone() const
        {
            return autoPtr<StochasticCollisionModel<CloudType>>
            (
                new TrajectoryCollision<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~TrajectoryCollision();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/StochasticCollision/TrajectoryCollision/TrajectoryCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
