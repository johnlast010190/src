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
    (c) 2012-2014 OpenFOAM Foundation

Class
    Foam::patchSeedSet

Description
    Initialises points on or just off patch

SourceFiles
    patchSeedSet.C

\*---------------------------------------------------------------------------*/

#ifndef patchSeedSet_H
#define patchSeedSet_H

#include "sampledSet/sampledSet/sampledSet.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "primitives/random/Random/Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class patchSeedSet Declaration
\*---------------------------------------------------------------------------*/

class patchSeedSet
:
    public sampledSet
{
    // Private data

        //- Patches to sample
        const labelHashSet patchSet_;

        //- Maximum number of patch faces to seed (if in random subset mode)
        const label maxPoints_;

        //- Random number generator (if maxPoints < num patch faces)
        autoPtr<Random> rndGenPtr_;

        //- Patch faces to seed selected based on nearness to supplied points
        const pointField selectedLocations_;


    // Private Member Functions

        //- Samples all points in sampleCoords.
        void calcSamples
        (
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        );

        //- Uses calcSamples to obtain samples. Copies them into *this.
        void genSamples();


public:

    //- Runtime type information
    TypeName("patchSeed");


    // Constructors

        //- Construct from dictionary
        patchSeedSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~patchSeedSet();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
