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
    (c) 2023 Esi Ltd.

Class
    Foam::foamMeshGenerator

Description
    Encapsulates generation of a foam mesh

SourceFiles
    foamMeshGenerator.C

\*---------------------------------------------------------------------------*/

#include "fvMeshTools/fvMeshTools.H"
#include "snappyHexMeshDriver/snappyRefineDriver.H"
#include "snappyHexMeshDriver/snappySnapDriver.H"
#include "snappyHexMeshDriver/snappyLayerDriver.H"

#ifndef foamMeshGenerator_H
#define foamMeshGenerator_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class foamMeshGenerator Declaration
\*---------------------------------------------------------------------------*/

class foamMeshGenerator
{

public:


private:

    // Private data

        //- reference to Time
        Time& runTime_;

    // Member Functions

        //- Generate optimisation setup for final stages
        dictionary finalOptimizeSetup
        (
            const dictionary& meshDict
        );

        //- Perform cfMesh opyimisation stage
        void cfOptimize
        (
            fvMesh& mesh,
            const dictionary& meshOptimDict
        );

        //- Extract triangulated surface from final meshed surface
        void extractTriSurface
        (
            const polyMesh& mesh,
            const labelHashSet& includePatches,
            const fileName& outFileName
        );

        //- Calculate distance used when merging points
        scalar getMergeDistance
        (
            const polyMesh& mesh,
            const boundBox& meshBb,
            const scalar mergeTol
        );

        //- Write out final mesh
        void writeMesh
        (
            const string& msg,
            meshRefinement& meshRefiner,

            const meshRefinement::debugType debugLevel,
            const meshRefinement::writeType writeLevel,

            bool writeSurf  = true,
            bool removeZeroSizedPatches  = true
        );

        //- Create tri surface from blockMesh patches (extrude/dual method)
        void blockMeshTriangulate
        (
            const polyMesh& mesh,
            dictionary& geometryDict,
            dictionary& refineDict,
            const bool baseCheck
        );

        //- Add VDB patches to the mesh
        void addVDBpatches
        (
            const dictionary dict,
            fvMesh& mesh
        );

        //- Create helper dictionary for re-snapping to
        // zoned internal boundaries during layer projection (dual/extrude)
        dictionary calcReSnapNonBoundaryZoneDict
        (
            const layerParameters& layerParams,
            const refinementSurfaces& surfaces,
            const dictionary& geometryDict
        );

        //- Create helper dictionary for re-snapping to
        // grown-up boundaries during layer projection (dual/extrude)
        dictionary calcReSnapGeometryDict
        (
            const layerParameters& layerParams,
            const refinementSurfaces& surfaces,
            const labelList& globalToMasterPatch,
            const labelList& globalToSlavePatch,
            const dictionary& geometryDict
        );

        //- Automatically create hierarchical coefficients
        void autoCreateDecomposition
        (
            const label nProcs,
            dictionary& decomposeDict
        );

        //- For calculating protected cell field for dynamic refinement
        void addProtectedCellsExtrude
        (
            const fvMesh& mesh,
            volScalarField& protectedCells
        );

        //- Create a starting blockMesh from blockMeshDict and re-partition
        // if running in parallel
        void runBlockMesh
        (
            const dictionary& meshDict,
            const word& regionName,
            const word& regionDir,
            autoPtr<fvMesh>& meshPtr,
            decompositionMethod& decomposer
        );

        //- Optional addition of layers
        bool addLayers
        (
            const dictionary& layerDict,
            const dictionary& motionDict,
            const meshControl& controller,
            const labelList& globalToMasterPatch,
            const labelList& globalToSlavePatch,
            const bool& overwrite,
            meshRefinement& meshRefiner,
            autoPtr<searchableSurfaces>& allGeometryPtr,
            hexReport& stats,
            decompositionMethod& decomposer,
            fvMeshDistribute& distributor,
            const label multiLayerIter
        );

public:

    // Constructors

        //- Construct from Time
        foamMeshGenerator(Time&);

    // Member Functions

        // Generate and return a mesh object
        autoPtr<fvMesh> generateMesh
        (
           const word& regionName,
           const word& regionDir,
           const fileName& dictPath,
           const fileName& decompDictFile,
           const bool& baseCheck,
           const bool& overwrite,
           const bool& writeIntermediate,
           const bool& writeConstant,
           const scalar& writeTime,
           const bool& checkGeometry,
           bool writeDict,
           bool meshWrite = true
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
