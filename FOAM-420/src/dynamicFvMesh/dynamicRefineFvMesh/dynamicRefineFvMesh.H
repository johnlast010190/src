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
    (c) 2011-2014 OpenFOAM Foundation
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.
    (c) 2016 Esi Ltd.

Class
    Foam::dynamicRefineFvMesh

Description
    A fvMesh with built-in refinement.

    Determines which cells to refine/unrefine and does all in update().

SourceFiles
    dynamicRefineFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicRefineFvMesh_H
#define dynamicRefineFvMesh_H

#include "dynamicFvMesh/dynamicFvMesh.H"
#include "polyTopoChange/hexRef8/hexRef8.H"
#include "containers/Lists/PackedList/PackedBoolList.H"
#include "primitives/bools/Switch/Switch.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "fvMesh/fvPatches/constraint/cyclicAMI/cyclicAMIFvPatch.H"
#include "decompositionMethod/decompositionMethod.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "regionSplit/regionSplit.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "geomDecomp/geomDecomp.H"
#include "fvMeshDistribute/fvMeshDistribute.H"
#include "turbulenceModel.H"
#include "fvMesh/wallDist/nearWallDist/nearWallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicRefineFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicRefineFvMesh
:
    public dynamicFvMesh
{
protected:

        //- Mesh cutting engine
        hexRef8 meshCutter_;

        //- Dump cellLevel for postprocessing
        Switch dumpLevel_;

        //- Fluxes to map
        HashTable<word> correctFluxes_;

        //- SurfaceFiels to map
        HashTable<word> mapSurfaceFields_;

        //- Number of refinement/unrefinement steps done so far.
        label nRefinementIterations_;

        //- Protected cells (usually since not hexes)
        PackedBoolList protectedCell_;

        //- protect all wall cells
        //  reduces issues with wallFunctions during dynamic refine
        Switch protectLayers_;

        //- switch to disable y correction in turbulence models
        Switch correctY_;

        //- protect cells near specicif patches
        wordList protectedPatches_;

        //- protect cells in specicif cellZone
        wordList protectedCellZones_;

        //- global number of protected cells
        label gnProtected_;

        //- global number of protected wall-cells
        label gnWallCellsProtected_;

        //- global number of protected patch-cells
        label gnPatchCellsProtected_;

        //- global number of protected cellZone-cells
        label gnCellZoneCellsProtected_;

        //- start refining at this time
        mutable scalar refineStartTime_;

        //- stop refining at this time
        mutable scalar refineStopTime_;

        //- allow refinement of dangling cells
        Switch danglingRefine_;

        //- if available read and maintain field of protected cells
        autoPtr<volScalarField> protectedCells_;

        //- refine mesh outer correctors
        bool refineMeshOuterCorr_;

        //- is first iteration
        bool isFirstIter_;

    // Protected Member Functions

        //- Count set/unset elements in packedlist.
        static label count(const PackedBoolList&, const unsigned int);

        //- Calculate cells that cannot be refined since would trigger
        //  refinement of protectedCell_ (since 2:1 refinement cascade)
        void calculateProtectedCells(PackedBoolList& unrefineableCell) const;

        //- find the protected cells since redistribution will
        //  change the local proc protectedCell_ list
        void findProtectedCells();

        //- Read the projection parameters from dictionary
        void readDict();

        //- Refine cells. Update mesh and fields.
        autoPtr<mapPolyMesh> refine(const labelList&);

        //- Unrefine cells. Gets passed in centre points of cells to combine.
        autoPtr<mapPolyMesh> unrefine(const labelList&);


        // Selection of cells to un/refine

            //- Calculates approximate value for refinement level so
            //  we don't go above maxCell
            scalar getRefineLevel
            (
                const label maxCells,
                const label maxRefinement,
                const scalar refineLevel,
                const scalarField&
            ) const;

            //- Get per cell max of connected point
            scalarField maxPointField(const scalarField&) const;

            //- Get point min of connected cell
            virtual scalarField minCellField(const volScalarField&) const;

            scalarField cellToPoint(const scalarField& vFld) const;

            scalarField error
            (
                const scalarField& fld,
                const scalar minLevel,
                const scalar maxLevel
            ) const;

            //- Select candidate cells for refinement
            virtual void selectRefineCandidates
            (
                const scalar lowerRefineLevel,
                const scalar upperRefineLevel,
                const scalarField& vFld,
                PackedBoolList& candidateCell
            ) const;

            //- Subset candidate cells for refinement
            virtual labelList selectRefineCells
            (
                const label maxCells,
                const label maxRefinement,
                const PackedBoolList& candidateCell
            ) const;

            //- Select points that can be unrefined.
            virtual labelList selectUnrefinePoints
            (
                const scalar unrefineLevel,
                const PackedBoolList& markedCell,
                const scalarField& pFld
            ) const;

            //- Extend markedCell with cell-face-cell.
            void extendMarkedCells(PackedBoolList& markedCell) const;

            //- Check all cells have 8 anchor points
            void checkEightAnchorPoints
            (
                PackedBoolList& protectedCell,
                label& nProtected
            ) const;

            //- recursive refinement history search
            label parentID(label p);

            //- only unrefine cells that have been refined
            void checkRefHistory
            (
                PackedBoolList& hasRefHistory
            ) const;

            bool walkNeighbours
            (
                const label& patchMaster,
                const label& patchSlave,
                labelList& agglom,
                point& avgCC,
                label& nAMICells
            );

            //- Refine cells with almost all sides refined
            //  alter marked cells analogous to snappeRefineDriver
            void danglingCellRefine
            (
                labelList& cellsToRefine,
                const label nFaces
            );

            //- interpolate to newly generated internal faces of a parent cell
            template <class T>
            void mapNewInternalFaces(const labelList& faceMap);

private:

        //- Disallow default bitwise copy construct
        dynamicRefineFvMesh(const dynamicRefineFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicRefineFvMesh&);

public:

    //- Runtime type information
    TypeName("dynamicRefineFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicRefineFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicRefineFvMesh();


    // Member Functions

        //- Direct access to the refinement engine
        const hexRef8& meshCutter() const
        {
            return meshCutter_;
        }

        //- Cells which should not be refined/unrefined
        const PackedBoolList& protectedCell() const
        {
            return protectedCell_;
        }

        //- Cells which should not be refined/unrefined
        PackedBoolList& protectedCell()
        {
            return protectedCell_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Template to update all volField boundaries
        template<class Type> void correctBoundaries();

    // Writing

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool valid
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
