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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2015-2018 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::polyPatch

Description
    polyPatch is the base class of:
    -the direct polyPatches (created by primitivePatch)
    -the indirect polyPatches (created by indirectPrimitivePatch)

    A patch is a list of labels that address the faces in the global face list.

    The patch can calculate its own edges based on the global faces.
    Patch also contains all addressing between the faces.

SourceFiles
    polyPatch.C
    polyPatchNew.C

\*---------------------------------------------------------------------------*/

#ifndef polyPatch_H
#define polyPatch_H

#include "meshes/Identifiers/patch/patchIdentifier.H"
#include "meshes/primitiveMesh/primitivePatch/primitivePatch.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declaration of friend functions and operators

class polyBoundaryMesh;
class polyPatch;
class polyTopoChange;
class PstreamBuffers;
class face;

Ostream& operator<<(Ostream&, const polyPatch&);

/*---------------------------------------------------------------------------*\
                           Class polyPatch Declaration
\*---------------------------------------------------------------------------*/

class polyPatch
:
    public patchIdentifier
{
    // Private data

protected:

        //- Start label of this patch in the polyMesh face list
        label start_;

    // Protected Member Functions

        // The polyPatch geometry initialisation is called by polyBoundaryMesh
        friend class polyBoundaryMesh;

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&) = 0;

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&) = 0;

        //- Correct patches after moving points
        virtual void movePoints(const pointField& p) = 0;

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&) = 0;

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField& p) = 0;

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&) = 0;

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&) = 0;

        virtual void updateGIB(){}

        //- Clear geometry
        virtual void clearGeom() = 0;

        //- Reset the patch name
        virtual void rename(const wordList& newNames) = 0;

        //- Reset the patch index
        virtual void reorder(const labelUList& oldToNewIndex) = 0;


public:

    //- Runtime type information
    TypeName("basePatch");

    //- Debug switch to disallow the use of genericPolyPatch
    static int disallowGenericPolyPatch;


    // Declare run-time constructor selection tables

        //for direct polypatches
        declareRunTimeSelectionTable
        (
            autoPtr,
            polyPatch,
            word,
            (
                const word& name,
                const label size,
                const label start,
                const label index,
                const polyBoundaryMesh& bm,
                const word& patchType
            ),
            (name, size, start, index, bm, patchType)
        );

        //for indirect polypatches
        declareRunTimeSelectionTable
        (
            autoPtr,
            polyPatch,
            indirect,
            (
                const word& name,
                const label size,
                const label zoneI,
                const word zPType,
                const label index,
                const polyBoundaryMesh& bm,
                const word& patchType
            ),
            (name, size, zoneI, zPType, index, bm, patchType)
        );

        //from dictionary
        declareRunTimeSelectionTable
        (
            autoPtr,
            polyPatch,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const label index,
                const polyBoundaryMesh& bm,
                const word& patchType
            ),
            (name, dict, index, bm, patchType)
        );


    // Constructors

        //- Construct from components directPatch
        polyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from components directPatch
        polyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& physicalType,
            const wordList& inGroups
        );

        //- Construct from components indirectPatch
        polyPatch
        (
            const word&,
            const label,
            const label start,
            const label,
            const word,
            const label,
            const polyBoundaryMesh&,
            const word&
        );

        //- Construct from dictionary
        polyPatch
        (
            const word& name,
            const label start,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        polyPatch(const polyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        polyPatch
        (
            const polyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct given the original patch and a map
        polyPatch
        (
            const polyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct as copy
        polyPatch(const polyPatch&);

        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const = 0;

    // Selectors

        //- Return a pointer to a new patch created on freestore from
        //  components - direct pp
        static autoPtr<polyPatch> New
        (
            const word& patchType,
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );


        //- Return a pointer to a new patch created on freestore from
        //  components - indirect pp
        static autoPtr<polyPatch> New
        (
            const word& patchType,
            const word& name,
            const label size,
            const label zoneId,
            const word  zPType,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Return a pointer to a new patch created on freestore from
        //  dictionary
        static autoPtr<polyPatch> New
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Return a pointer to a new patch created on freestore from
        //  dictionary
        static autoPtr<polyPatch> New
        (
            const word& patchType,
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm
        );


    //- Destructor


    // Member Functions

        //- Return boundaryMesh reference
        virtual const polyBoundaryMesh& boundaryMesh() const = 0;

        //- Return true if this patch is geometrically coupled (i.e. faces and
        //  points correspondence)
        virtual bool coupled() const
        {
            return false;
        }

        //- Return start label of this patch in the polyMesh face list
        label start() const {
            return this->start_;
        }

        label& start() {
            return this->start_;
        }

        virtual label localToGlobal(label&) const = 0;

        //- Return true if the given type is a constraint type
        static bool constraintType(const word& pt);

        //- Return a list of all the constraint patch types
        static wordList constraintTypes();

        //- Extract face cell data
        template<class T>
        const UIndirectList<T> patchInternalList
        (
            const UList<T>& internalValues
        ) const
        {
            return UIndirectList<T>(internalValues, faceCells());
        }

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const = 0;


        //- remove the contributions from the matrix
        //- used in indirect p
        virtual void GIBAgglomerate
        (
            scalarField&
        ) const=0;

        //- virtual functions
        virtual const List<scalar>::subList patchSlice
        (
            const Field<scalar>&
        ) const =0;
        virtual const List<vector>::subList patchSlice
        (
            const Field<vector>&
        ) const =0;
        virtual const List<scalar>::subList patchSlice
        (
            const List<scalar>&
        ) const =0;
        virtual const List<vector>::subList patchSlice
        (
            const List<vector>&
        ) const =0;
        virtual const List<vector> indirectPatchSlice
        (
            const Field<vector>&
        ) const =0;
        virtual const List<scalar> indirectPatchSlice
        (
            const Field<scalar>&
        ) const =0;

        // Geometric data; point list required

            //- Return face centres
            virtual const vectorField::subField faceCentres() const = 0;

            //- Return face areas
            virtual const vectorField::subField faceAreas() const = 0;

            //- Return face area magnitudes
            virtual const scalarField::subField magFaceAreas() const = 0;

            //- Return face neighbour cell centres
            virtual tmp<vectorField> faceCellCentres() const = 0;


        // Addressing into mesh

            //- Return face-cell addressing
            virtual const labelUList& faceCells() const = 0;

            //- Return global edge index for local edges
            virtual const labelList& meshEdges() const = 0;

            //- Clear addressing
            virtual void clearAddressing() = 0;

        // Pure functions from primitive mesh
            virtual label size() const=0;
            virtual const labelList& boundaryPoints() const=0;
            virtual const labelListList& faceFaces() const=0;
            virtual const labelListList& edgeFaces() const=0;
            virtual const labelListList& faceEdges() const=0;
            virtual const labelListList& pointEdges() const=0;
            virtual const labelListList& pointFaces() const=0;
            virtual label nPoints() const=0;
            virtual label nEdges() const=0;
            virtual const edgeList& edges() const=0;
            virtual label nInternalEdges() const=0;
            virtual bool isInternalEdge(const label edgeI) const=0;
            virtual const List<face>& localFaces() const=0;
            virtual const labelList& meshPoints() const=0;
            virtual const Map<label>& meshPointMap() const=0;
            virtual const Field<point>& localPoints() const=0;
            virtual const Field<point>& faceNormals() const=0;
            virtual const Field<point>& pointNormals() const = 0;
            virtual const Field<point>& points() const = 0;
            virtual bool empty() const = 0;
            virtual const labelListList& edgeLoops() const = 0;
            virtual const Field<scalar>& faceCurvature() const = 0;
            virtual const Field<scalar>& pointCurvature() const = 0;

        // Other patch operations
            virtual label whichFace(const label) const = 0;
            virtual label whichPoint(const label gp) const = 0;
            virtual label whichEdge(const edge&) const = 0;


        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder(PstreamBuffers&, const primitivePatch&) const = 0;

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
        ) const = 0;

        //- For dynamic mesh cases - return true if this patch will change the
        //- topology
        virtual bool changeTopology() const
        {
            return false;
        }

        //- Collect topology changes in a polyTopoChange object
        virtual bool setTopology(polyTopoChange&)
        {
            return false;
        }


    // Member operators
        virtual void operator=(const polyPatch& p) = 0;
        virtual const face& operator[](const label) const = 0;
        virtual face& operator[](const label) = 0;

    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const polyPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
