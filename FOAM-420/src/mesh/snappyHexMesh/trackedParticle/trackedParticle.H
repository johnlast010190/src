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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2013 OpenFOAM Foundation

Class
    Foam::trackedParticle

Description
    Particle class that marks cells it passes through. Used to mark cells
    visited by feature edges.

SourceFiles
    trackedParticle.C

\*---------------------------------------------------------------------------*/

#ifndef trackedParticle_H
#define trackedParticle_H

#include "particle/particle.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class trackedParticleCloud;


// Forward declaration of friend functions and operators

class trackedParticle;

Ostream& operator<<(Ostream&, const trackedParticle&);


/*---------------------------------------------------------------------------*\
                     Class trackedParticle Declaration
\*---------------------------------------------------------------------------*/

class trackedParticle
:
    public particle
{
    // Private data

        //- end point to track to
        point end_;

        //- level of this particle
        label level_;

        //- passive label
        label i_;

        //- passive label
        label j_;

        //- particle on external patch
        bool onPatch_;

        //- particle on ftr line
        bool onFtr_;

public:

    friend class Cloud<trackedParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<Cloud<trackedParticle>>
    {
        labelList& maxLevel_;

    public:


        // Constructors

            trackingData(Cloud<trackedParticle>& cloud, labelList& maxLevel)
            :
                particle::TrackingData<Cloud<trackedParticle>>(cloud),
                maxLevel_(maxLevel)
            {}


        // Member functions

            labelList& maxLevel()
            {
                return maxLevel_;
            }
    };



    // Constructors

        //- Construct from components
        trackedParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label tetFacei,
            const label tetPtI,
            const point& end,
            const label level,
            const label i,
            const label j
        );

        //- Construct from Istream
        trackedParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new trackedParticle(*this));
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

            autoPtr<trackedParticle> operator()(Istream& is) const
            {
                return autoPtr<trackedParticle>
                (
                    new trackedParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        //- point to track to
        point& end()
        {
            return end_;
        }

        //- transported label
        label i() const
        {
            return i_;
        }

        //- transported label
        label& i()
        {
            return i_;
        }

        //- transported label
        label& j()
        {
            return j_;
        }

        //- particle stuck to patch
        bool& onPatch()
        {
            return onPatch_;
        }

        //- particle on feature line
        bool& onFtr()
        {
            return onFtr_;
        }

        // Tracking

            //- Track all particles to their end point
            bool move(trackingData&, const scalar);


            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch&,
                trackingData& td,
                const label patchi,
                const scalar trackFraction,
                const tetIndices& tetIs
            );

            //- Overridable function to handle the particle hitting a wedge
            void hitWedgePatch
            (
                const wedgePolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //  symmetry plane
            void hitSymmetryPlanePatch
            (
                const symmetryPlanePolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //  symmetry patch
            void hitSymmetryPatch
            (
                const symmetryPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a cyclic
            void hitCyclicPatch
            (
                const cyclicPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //- processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                const wallPolyPatch&,
                trackingData& td,
                const tetIndices&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                trackingData& td
            );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const trackedParticle&);
};


template<>
constexpr bool contiguous<trackedParticle>()
{
    return true;
}

//template<>
//void Cloud<trackedParticle>::readFields();
//
//template<>
//void Cloud<trackedParticle>::writeFields() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
