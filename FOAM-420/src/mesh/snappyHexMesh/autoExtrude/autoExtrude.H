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
    (c) 2019 Esi Ltd.

Application
    Foam::autoExtrude

Description
    Add layers of cells to all boundary faces

\*---------------------------------------------------------------------------*/

#ifndef autoExtrude_H
#define autoExtrude_H

#include "meshRefinement/meshRefinement.H"
#include "regionSplit2D/regionSplit2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class mapPolyMesh;
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class autoExtrude Declaration
\*---------------------------------------------------------------------------*/

class autoExtrude
{
public:

    // Public data

        //- extrusion type enumeration
        enum extType
        {
            extDirection,     //Fixed direction
            extNormal,        //Patch normal direction
            extAverage,       //Average patch normal direction
            extTarget         //Target patch mapping
        };

        static const NamedEnum<extType, 4> extTypeNames;

        //- extrusion type enumeration
        enum extAction
        {
            ADD,     //Add to mesh
            NEW,        //Create new mesh
            SPLIT       //Split into new region
        };

        static const NamedEnum<extAction, 3> extActionNames;

private:

    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- Reference to decomposition method
        decompositionMethod& decomposer_;

        //- Reference to mesh distribution engine
        fvMeshDistribute& distributor_;

    // Private Member Functions

public:

    //- Runtime type information
    ClassName("autoExtrude");

    // Constructors

        //- Construct from mesh.
        autoExtrude
        (
            meshRefinement& meshRefiner,
            decompositionMethod& decomposer,
            fvMeshDistribute& distributor
        );

    // Member Functions

        // Helper

            void addNewProcessorPatches
            (
                const indirectPrimitivePatch& pp,
                labelList& edgePatchID,
                labelList& edgeZoneID,
                boolList& edgeFlip,
                labelList& inflateFaceID
            );

            void calcRegionsProperties
            (
                const indirectPrimitivePatch& pp,
                const regionSplit2D& regionIDs,

                scalarField& extrudeFCH,
                pointField& regionCentroids,
                pointField& regionNormals,
                scalarField& regionAreas
            );

        // Access

        // extrude boundary faces

            //Filter faces that should not be extruded
            void filterExtrudePatch
            (
                const labelList& patchIDs
            );

            // Calculate based on input displacement normal
            bool calculateDirectionDisplacement
            (
                const dictionary& edict,
                const indirectPrimitivePatch& extrudePatch,
                const regionSplit2D& sourceRegionIDs,
                pointField& displacement
            );

            // Calculate based on calculated displacement normal
            bool calculateNormalDisplacement
            (
                const dictionary& edict,
                const indirectPrimitivePatch& extrudePatch,
                const regionSplit2D& sourceRegionIDs,
                pointField& displacement
            );

            // Calculate based on calculated average displacement normal
            bool calculateAverageNormalDisplacement
            (
                const dictionary& edict,
                const indirectPrimitivePatch& extrudePatch,
                const regionSplit2D& sourceRegionIDs,
                pointField& displacement
            );

            scalar calculateStretch
            (
                const scalar lheight,
                const scalar fch,
                const label nLayers
            );

            // Calculate based on input target patch (planar and translational)
            bool calculateTargetDisplacement
            (
                const dictionary& edict,
                const indirectPrimitivePatch& extrudePatch,
                const labelHashSet& sourceSet,
                const regionSplit2D& sourceRegionIDs,
                const label nLayers,
                pointField& displacement,
                scalarField& expansionRatio
            );

            // Add optional zoning for final patch extrusion
            void addZones
            (
                const dictionary& edict,
                labelList& cellZoneSource,
                List<labelPair>& faceZoneSource
            );

            // Before extruding flatten source (and target) patches
            void flattenPatch
            (
                const indirectPrimitivePatch& pp
            );

            //Clean-up free-standing zone faces
            void cleanupFaceZones();

            void extrudeSelected
            (
                const dictionary& extrudeDict
            );

            void extrudeZeroSizedLayer
            (
                const refinementParameters& refineParams
            );

            void removeUnusedPoints();

            void optimize
            (
                const refinementParameters& refineParams
            );

            void preBalance
            (
                const dictionary& extrudeDict,
                const bool updateSurf
            );

            void balance
            (
                const scalar snapWeights,
                const bool updateSurf
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
