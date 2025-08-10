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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2015-2018 OpenCFD Ltd.
    (c) 2015-2022 Esi Ltd.

Class
    Foam::directPolyPatch

Description
    A direct polyPatch is created by PrimitivePatch (standard behaviour)
    A patch is a list of labels that address the faces in the global face list.

    The patch can calculate its own edges based on the global faces.
    Patch also contains all addressing between the faces.

SourceFiles
    directPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef directPolyPatch_H
#define directPolyPatch_H

#include "meshes/primitiveMesh/primitivePatch/primitivePatch.H"
#include "meshes/polyMesh/polyPatches/polyPatch/polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declaration of friend functions and operators

class directPolyPatch;
Ostream& operator<<(Ostream&, const directPolyPatch&);

/*---------------------------------------------------------------------------*\
                           Class directPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class directPolyPatch
:
    public polyPatch,
    public primitivePatch
{
    // Private data

        //- Reference to boundary mesh
        const polyBoundaryMesh& boundaryMesh_;


        // Demand-driven private data

            //- face-cell addressing
            mutable labelList::subList* faceCellsPtr_;

            //- Global edge addressing
            mutable labelList* mePtr_;


protected:

    // Protected Member Functions

        // The directPolyPatch geometry initialisation is called by polyBoundaryMesh
        friend class polyBoundaryMesh;

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&)
        {}

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&)
        {}

        //- Correct patches after moving points
        virtual void movePoints(const pointField& p);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&)
        {}

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField& p);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&)
        {}

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);

        virtual void updateGIB();

        //- Clear geometry
        virtual void clearGeom();

        //- Reset the patch name
        virtual void rename(const wordList& newNames);

        //- Reset the patch index
        virtual void reorder(const labelUList& oldToNewIndex);


public:

    //- Runtime type information
    TypeName("patch");

    // Constructors

        //- Construct from components
        directPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from components
        directPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& physicalType,
            const wordList& inGroups
        );

        //- Construct from dictionary
        directPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        directPolyPatch(const directPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        directPolyPatch
        (
            const directPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        directPolyPatch
        (
            const polyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct given the original patch and a map
        directPolyPatch
        (
            const directPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct as copy
        directPolyPatch(const directPolyPatch&);

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new directPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<directPolyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<directPolyPatch>
            (
                new directPolyPatch(*this, bm, index, newSize, newStart)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<directPolyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<directPolyPatch>
            (
                new directPolyPatch(*this, bm, index, mapAddressing, newStart)
            );
        }


    //- Destructor
    virtual ~directPolyPatch();


    // Member Functions

        //- The offset where this patch starts in the boundary face list
        //  The value is the same as  patch.start() - mesh.nInternalFaces()
        label offset() const;


        virtual label localToGlobal(label& fi) const
        {
            return this->start_+fi;
        }

        //- Return boundaryMesh reference
        const polyBoundaryMesh& boundaryMesh() const;

        //- Return true if this patch is geometrically coupled (i.e. faces and
        //  points correspondence)
        virtual bool coupled() const
        {
            return false;
        }

        //- Return true if the given type is a constraint type
        static bool constraintType(const word& pt);

        //- Return a list of all the constraint patch types
        static wordList constraintTypes();

        virtual void GIBAgglomerate
        (
            scalarField&
        ) const
        {}

        //- Slice list to patch
        template<class T>
        const typename List<T>::subList patchSlice(const UList<T>& l) const
        {
            return typename List<T>::subList(l, this->size(), start_);
        }

        //- Slice Field to patch
        template<class T>
        const typename Field<T>::subField patchSlice(const Field<T>& l) const
        {
            return typename Field<T>::subField(l, this->size(), start_);
        }

        virtual const List<scalar>::subList patchSlice
        (
            const Field<scalar>& l
        ) const
        {
            return List<scalar>::subList(l, this->size(), start_);
        }
        virtual const List<vector>::subList patchSlice
        (
            const Field<vector>& l
        ) const
        {
            return List<vector>::subList(l, this->size(), start_);
        }
        virtual const List<scalar>::subList patchSlice
        (
            const List<scalar>& l
        ) const
        {
            return List<scalar>::subList(l, this->size(), start_);
        }
        virtual const List<vector>::subList patchSlice
        (
            const List<vector>& l
        ) const
        {
            return List<vector>::subList(l, this->size(), start_);
        }
        virtual const List<vector> indirectPatchSlice
        (
            const Field<vector>& l
        ) const
        {
            NotImplemented;
        }
        virtual const List<scalar> indirectPatchSlice
        (
            const Field<scalar>& l
        ) const
        {
            NotImplemented;
        }

        //- Write the directPolyPatch data as a dictionary
        virtual void write(Ostream&) const;


        // Geometric data; point list required

            //- Return face centres
            virtual const vectorField::subField faceCentres() const;

            //- Return face areas
            virtual const vectorField::subField faceAreas() const;

            //- Return face area magnitudes
            virtual const scalarField::subField magFaceAreas() const;

            //- Return face neighbour cell centres
            virtual tmp<vectorField> faceCellCentres() const;


        // Addressing into mesh

            //- Return face-cell addressing
            virtual const labelUList& faceCells() const;

            //- Return global edge index for local edges
            virtual const labelList& meshEdges() const;

            //- Clear addressing
            virtual void clearAddressing();

        // Pure functions from primitivePatch
            virtual label size() const
            {
                return this->primitivePatch::size();
            }
            virtual const labelList& boundaryPoints() const
            {
                return this->primitivePatch::boundaryPoints();
            }
            virtual const labelListList& faceFaces() const
            {
                return this->primitivePatch::faceFaces();
            }

            virtual const labelListList& edgeFaces() const
            {
                return this->primitivePatch::edgeFaces();
            }
            virtual const labelListList& faceEdges() const
            {
                return this->primitivePatch::faceEdges();
            }
            virtual const labelListList& pointEdges() const
            {
                return this->primitivePatch::pointEdges();
            }
            virtual const labelListList& pointFaces() const
            {
                return this->primitivePatch::pointFaces();
            }
            virtual label nPoints() const
            {
                return this->primitivePatch::nPoints();
            }
            virtual label nEdges() const
            {
                return this->primitivePatch::nEdges();
            }
            virtual const edgeList& edges() const
            {
                return this->primitivePatch::edges();
            }
            virtual label nInternalEdges() const
            {
                return this->primitivePatch::nInternalEdges();
            }
            virtual bool isInternalEdge(const label edgeI) const
            {
                return this->primitivePatch::isInternalEdge(edgeI);
            }
            virtual const List<face>& localFaces() const
            {
                return this->primitivePatch::localFaces();
            }
            virtual const labelList& meshPoints() const
            {
                return this->primitivePatch::meshPoints();
            }
            virtual const Map<label>& meshPointMap() const
            {
                return this->primitivePatch::meshPointMap();
            }
            virtual const Field<point>& localPoints() const
            {
                return this->primitivePatch::localPoints();
            }
            const Field<point>& faceNormals() const
            {
                return this->primitivePatch::faceNormals();
            }
            const Field<point>& pointNormals() const
            {
                return this->primitivePatch::pointNormals();
            }
            const Field<point>& points() const
            {
                return this->primitivePatch::points();
            }
            virtual bool empty() const
            {
                return this->primitivePatch::empty();
            }
            virtual const labelListList& edgeLoops() const
            {
                return this->primitivePatch::edgeLoops();
            }
            virtual const Field<scalar>& faceCurvature() const
            {
                return this->primitivePatch::faceCurvature();
            }
            virtual const Field<scalar>& pointCurvature() const
            {
                return this->primitivePatch::pointCurvature
                (
                    boolList(this->size(), false)
                );
            }

        // Other patch operations

            //- Return label of face in patch from global face label
            virtual label whichFace(const label l) const
            {
                return l - start_;
            }

            virtual label whichPoint(const label gp) const
            {
                return this->primitivePatch::whichPoint(gp);
            }

            virtual label whichEdge(const edge& ed) const
            {
                return this->primitivePatch::whichEdge(ed);
            }


        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder(PstreamBuffers&, const primitivePatch&) const{}

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            PstreamBuffers&,
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const
        {
            return false;
        }


    // Member operators

        //- Assignment
        virtual void operator=(const polyPatch&);
        virtual const face& operator[](const label) const;
        virtual face& operator[](const label);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const directPolyPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
