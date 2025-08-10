/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM® : Professional Open-source CFD
|   o   O   o    |  Version : 4.2.0
|    o     o     |  Copyright © 2015 ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM® <http://www.openfoam.org/>.

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

Portions Copyright © 2011 OpenFOAM Foundation.

Class
    Foam::regionSplitPatch

Description
    This class separates the indirect primitive patch into distinct unconnected
    regions, each of which is then given a label according to globalNumbering().

SourceFiles
    regionSplitPatch.C

\*---------------------------------------------------------------------------*/

#ifndef regionSplitPatch_H
#define regionSplitPatch_H

#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "primitives/bools/lists/boolList.H"
#include "meshes/MeshObject/MeshObject.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                         Class regionSplitPatch Declaration
\*---------------------------------------------------------------------------*/

class regionSplitPatch
:
    public MeshObject<polyMesh, TopologicalMeshObject, regionSplitPatch>,
    public labelList
{
    // Private data


        mutable autoPtr<globalIndex> globalNumberingPtr_;

        const indirectPrimitivePatch& pp_;

        //- Given a seed cell label, fill cellRegion/faceRegion with markValue
        //  for contiguous region around it
        void fillSeedMask
        (
            labelList& faceRegion,
            labelList& edgeRegion,
            const label seedFaceID,
            const label markValue
        ) const;

        //- Calculate local region split. Return number of regions.
        label calcLocalRegionSplit
        (
            const boolList& blockedEdge,
            labelList& faceRegion
        ) const;

        //- Calculate global region split. Return globalIndex.
        autoPtr<globalIndex> calcRegionSplit
        (
            const bool doGlobalRegions,
            const boolList& blockedEdge,
            labelList& faceRegion
        ) const;

public:

    //- Runtime type information
    ClassName("regionSplitPatch");

    // Constructors

        //- Construct from mesh primitive patch and whether edge is blocked
        //  NOTE: blockedEdge has to be consistent across coupled faces!
        regionSplitPatch
        (
            const polyMesh&,
            const indirectPrimitivePatch&,
            const boolList& blockedEdge,
            const bool doGlobalRegions = Pstream::parRun()
        );

    // Member Functions

        //- Return global region numbering
        const globalIndex& globalNumbering() const
        {
            return globalNumberingPtr_();
        }

        //- Return local number of regions
        label nLocalRegions() const
        {
            return globalNumbering().localSize(Pstream::myProcNo());
        }

        //- Return total number of regions
        label nRegions() const
        {
            return globalNumbering().size();
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
