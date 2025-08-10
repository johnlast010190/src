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
    (c) 2019 Esi Ltd.

Class
    Foam::cuttingTriSurface

Description
    Constructs triSurface through mesh.

    No attempt at resolving degenerate cases. Since the cut faces are
    usually quite ugly, they will always be triangulated.

Note
    When the cutting triSurface coincides with a mesh face, the cell edge on the
    positive side of the triSurface is taken.

SourceFiles
    cuttingTriSurface.C

\*---------------------------------------------------------------------------*/

#ifndef cuttingTriSurface_H
#define cuttingTriSurface_H

#include "triSurface/triSurface.H"
#include "indexedOctree/treeDataTriSurface.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "meshes/meshShapes/face/faceList.H"
#include "MeshedSurface/MeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class primitiveMesh;

/*---------------------------------------------------------------------------*\
                       Class cuttingTriSurface Declaration
\*---------------------------------------------------------------------------*/

class cuttingTriSurface
:
    public MeshedSurface<face>
{
protected:

        //- Reference to surface to work on
        triSurface surface_;

private:

    //- Private typedef for convenience
    typedef MeshedSurface<face> MeshStorage;

    // Private data

        //- Optional tolerance to use in searches
        scalar tolerance_;

        //- List of cells cut by the triSurface
        labelList meshCells_;

        //- Octree for searches
        mutable autoPtr<indexedOctree<treeDataTriSurface>> treePtr_;

    // Private Member Functions

        //- Demand driven construction of the octree
        const indexedOctree<treeDataTriSurface>& tree() const;

        //- Determine mesh edges cutting surface
        void calcEdgeCuts
        (
            const primitiveMesh&,
            const bool planarCut,
            scalarField& cutEdges
        );

        //- Determine cut cells, possibly restricted to a list of cells
        void calcCutCells
        (
            const primitiveMesh&,
            const scalarField& cutEdges,
            const labelUList& cellIdLabels = labelUList::null()
        );

        //- Determine intersection points (cutPoints).
        void intersectEdges
        (
            const primitiveMesh&,
            const scalarField& cutEdges,
            List<label>& edgePoint
        );

        //- Walk circumference of cell, starting from startEdgeI crossing
        //  only cut edges. Record cutPoint labels in faceVerts.
        static bool walkCell
        (
            const primitiveMesh&,
            const labelUList& edgePoint,
            const label celli,
            label startEdgeI,
            DynamicList<label>& faceVerts
        );

        //- Determine cuts for all cut cells.
        void walkCellCuts
        (
            const primitiveMesh& mesh,
            const bool triangulate,
            const labelUList& edgePoint
        );


protected:

    // Constructors

        //- Construct triSurface description without cutting
        cuttingTriSurface(const triSurface&);


    // Protected Member Functions

        //- Recut mesh with existing triSurface, restricted to a list of cells
        void reCut
        (
            const primitiveMesh&,
            const bool triangulate,
            const bool planarCut = false,
            const labelUList& cellIdLabels = labelUList::null()
        );

        //- Remap action on triangulation or cleanup
        virtual void remapFaces(const labelUList& faceMap);


public:

    // Constructors

        //- Construct from triSurface and mesh reference,
        //  possibly restricted to a list of cells
        cuttingTriSurface
        (
            const triSurface&,
            const primitiveMesh&,
            const bool triangulate,
            const bool planarCut = false,
            const labelUList& cellIdLabels = labelUList::null()
        );


    // Member Functions

        //- Return reference to the surface.
        const triSurface& surface() const
        {
            return surface_;
        }

        //- Return List of cells cut by the triSurface
        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- Return true or false to question: have any cells been cut?
        bool cut() const
        {
            return meshCells_.size();
        }

        //- Sample the cell field
        template<class Type>
        tmp<Field<Type>> sample(const Field<Type>&) const;

        template<class Type>
        tmp<Field<Type>> sample(const tmp<Field<Type>>&) const;


    // Member Operators

        void operator=(const cuttingTriSurface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surface/cuttingTriSurface/cuttingTriSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
