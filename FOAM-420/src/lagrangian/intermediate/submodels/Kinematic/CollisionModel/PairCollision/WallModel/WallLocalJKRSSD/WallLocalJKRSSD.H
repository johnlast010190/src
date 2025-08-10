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
    (c) 2017 Esi Ltd.

Class
    Foam::WallLocalJKRSSD

Description
    JKR Ahesion Collision model from
    K. L. Johnson, K. Kendall, and A. D. Roberts
    "Surface Energy and Contact of Elastic Solids"
    proceedings of the Royal Society of London. Serias A.
    mathematical and Physical Sciences. Vol. 324, No 1558
    (Sept. 8, 1971), pp 301-313

    - Assumes normal contact of frictionless perticles (i.e. no tangential influences)
    - Approaches Hertz model as gamma approaches zero

\*---------------------------------------------------------------------------*/

#ifndef WallLocalJKRSSD_H
#define WallLocalJKRSSD_H

#include "submodels/Kinematic/CollisionModel/PairCollision/WallModel/WallModel/WallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                Class WallLocalJKRSSD Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class WallLocalJKRSSD
:
    public WallModel<CloudType>
{
    // Private data

        //- Effective Young's modulus value
        scalarList Estar_;

        //- Effective shear modulus value
        scalarList Gstar_;

        //- alpha-coefficient, related to coefficient of restitution
        scalarList alpha_;

        //- Coefficient of friction in for tangential sliding
        scalarList mu_;

        //- work of adhesion
        //- assumed identical for both components
        //- units of [kg/s^2]
        scalarList gamma_;

        //- Mapping the patch index to the model data
        labelList patchMap_;

        //- Index of the maximum value of Estar_
        label maxEstarIndex_;

        //- The number of steps over which to resolve the minimum
        //  harmonic approximation of the collision period
        scalar collisionResolutionSteps_;

        //- Volume factor for determining the equivalent size of a
        //  parcel where nParticles is not 1.  The equivalent size of
        //  the parcel is
        //      parcelEquivVolume = volumeFactor*nParticles*p.volume()
        //  so
        //      parcelEquivD = cbrt(volumeFactor*nParticles)*p.d()
        //  + When volumeFactor = 1, the particles are compressed
        //    together so that the equivalent volume of the parcel is
        //    the sum of the constituent particles
        //  + When volumeFactor = 3*sqrt(2)/pi, the particles are
        //    close packed, but uncompressed.
        //  + When volumeFactor > 3*sqrt(2)/pi, the particles loosely
        //    grouped.
        // 3*sqrt(2)/pi = 1.350474 is the volume factor for close
        // packing, i.e pi/(3*sqrt(2)) is the maximum close packing
        // factor
        scalar volumeFactor_;

        //- Switch to control use of equivalent size particles.  Used
        //  because the calculation can be very expensive.
        bool useEquivalentSize_;

        //- print out overlap and force infermation
        bool verbose_;

    // Private Member Functions

        //- Find the appropriate properties for determining the minimum
        //- allowable timestep
        void findMinMaxProperties
        (
            scalar& rMin,
            scalar& rhoMax,
            scalar& vMagMax
        ) const;

        //- Calculate the wall interaction for a parcel at a given site
        void evaluateWall
        (
            typename CloudType::parcelType& p,
            const point& site,
            const WallSiteData<vector>& data,
            scalar pREff,
            bool cohesion
        ) const;


public:

    //- Runtime type information
    TypeName("wallLocalJKRSSD");


    // Constructors

        //- Construct from dictionary
        WallLocalJKRSSD(const dictionary& dict, CloudType& cloud);


    //- Destructor
    virtual ~WallLocalJKRSSD();


    // Member Functions

        //- Return the volumeFactor
        inline scalar volumeFactor() const
        {
            return volumeFactor_;
        }

        //- Return the effective radius for a particle for the model
        virtual scalar pREff(const typename CloudType::parcelType& p) const;

        //- Whether the WallModel has a timestep limit that will
        //  require subCycling
        virtual bool controlsTimestep() const;

        //- For WallModels that control the timestep, calculate the
        //  number of subCycles needed to satisfy the minimum
        //  allowable timestep
        virtual label nSubCycles() const;

        //- Calculate the wall interaction for a parcel
        virtual void evaluateWall
        (
            typename CloudType::parcelType& p,
            const List<point>& flatSitePoints,
            const List<WallSiteData<vector>>& flatSiteData,
            const List<point>& sharpSitePoints,
            const List<WallSiteData<vector>>& sharpSiteData
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "submodels/Kinematic/CollisionModel/PairCollision/WallModel/WallLocalJKRSSD/WallLocalJKRSSD.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
