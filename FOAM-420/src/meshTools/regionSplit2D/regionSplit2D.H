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
    (c) 2015 OpenCFD Ltd.

Class
    Foam::regionSplit2D

Description
    Splits a patch into regions based on a mask field.  Result is a globally
    consistent label list of region index per patch face.

SourceFiles
    regionSplit2D.C

\*---------------------------------------------------------------------------*/

#ifndef regionSplit2D_H
#define regionSplit2D_H

#include "containers/Lists/DynamicList/DynamicList.H"
#include "primitives/bools/lists/boolList.H"
#include "primitives/ints/lists/labelList.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "algorithms/PatchEdgeFaceWave/patchEdgeFaceRegion.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faceZone;
class polyMesh;

/*---------------------------------------------------------------------------*\
                        Class regionSplit2D Declaration
\*---------------------------------------------------------------------------*/

class regionSplit2D
:
    public labelList
{
    // Private data

        //- Number of regions
        label nRegions_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionSplit2D(const regionSplit2D&);

        //- Disallow default bitwise assignment
        void operator=(const regionSplit2D&);


public:

    // Constructors

        //- Construct from mesh and list of blocked faces
        regionSplit2D
        (
            const polyMesh& mesh,
            const indirectPrimitivePatch& patch,
            const boolList& blockedFaces,
            const label offset = 0
        );


    //- Destructor
    ~regionSplit2D();


    // Member Functions

        //- Return the global number of regions
        label nRegions() const
        {
            return nRegions_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
