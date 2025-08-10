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
    (c) 2015-2024 Esi Ltd.

Class
    Foam::indirectPolyPatch

Description
    An indirect polyPatch is created by indirectPrimitivePatch based on a
    faceZone
    A patch is a list of labels that address the faces in the global face list.

    The patch can calculate its own edges based on the global faces.
    Patch also contains all addressing between the faces.

SourceFiles
    indirectPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef indirectPolyPatch_H
#define indirectPolyPatch_H

#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "meshes/polyMesh/polyPatches/polyPatch/polyPatch.H"
#include "fields/Fields/Field/SubField.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declaration of friend functions and operators

class indirectPolyPatch;
Ostream& operator<<(Ostream&, const indirectPolyPatch&);

/*---------------------------------------------------------------------------*\
                           Class indirectPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class indirectPolyPatch
:
    public polyPatch,
    public indirectPrimitivePatch
{

    // Private data

        //- Name of neighbour patch
        const word nbrPatchName_;

        //- Addressing from localToGlobal - based on the faceZone
        //  In boundary elements fAddr_!= faceZone
        mutable labelList fAddr_;

        //- flipmap
        mutable boolList fm_;

        //- Index of neighbour patch
        mutable label nbrPatchID_;

        label zoneId_;

        word indirectPolyPatchType_;

        //- Reference to boundary mesh
        const polyBoundaryMesh& boundaryMesh_;


        // Demand-driven private data

            //- face-cell addressing
            mutable labelList* faceCellsPtr_;

            //- Global edge addressing
            mutable labelList* mePtr_;

            mutable vectorField* faceCentresPtr_;

            mutable vectorField* faceAreasPtr_;

            mutable scalarField* magFaceAreasPtr_;

    // Private Member Functions

        //- create gib faceZone
        label createGIBFaceZone
        (
            const polyBoundaryMesh&,
            const word
        );

        //- modify the faceZone addressing for boundary faces
        labelList calculateAddressing
        (
            const word&,
            const polyBoundaryMesh&,
            const faceZone&
        );

        //- modify the faceZone addressing for boundary faces
        boolList calculateFm
        (
            const word&,
            const polyBoundaryMesh&,
            const faceZone&
        );


protected:

    // Protected Member Functions

        // The indirectPolyPatch geometry initialisation is called by polyBoundaryMesh
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
    TypeName("indirect");

    // Constructors

        //- Construct from components
        indirectPolyPatch
        (
            const word&,
            const label,
            const label,
            const word,
            const label,
            const polyBoundaryMesh&,
            const word&
        );

        //- Construct from dictionary
        indirectPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        indirectPolyPatch(const indirectPolyPatch&, const polyBoundaryMesh&);

        //- Construct as copy, resetting the faceZone
        indirectPolyPatch
        (
            const indirectPolyPatch&,
            const polyBoundaryMesh&,
            const label&
        );

        //- Construct as copy
        indirectPolyPatch(const indirectPolyPatch&);

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new indirectPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting thefaceZone
        //  and boundary mesh
        virtual autoPtr<indirectPolyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label& fzI
        ) const
        {
            return autoPtr<indirectPolyPatch>
            (
                new indirectPolyPatch(*this, bm, fzI)
            );
        }

    //- Destructor
    virtual ~indirectPolyPatch();


    // Member Functions

        const word& nbrPatchName() const
        {
            return nbrPatchName_;
        }

        //- Neighbour patchID.
        virtual label nbrPatchID() const;

        virtual bool owner() const
        {
            return index() < nbrPatchID();
        }

        const indirectPolyPatch& nbrPatch() const
        {
            const polyPatch& pp = this->boundaryMesh()[nbrPatchID()];
            return refCast<const indirectPolyPatch>(pp);
        }

        virtual label localToGlobal(label& fi) const
        {
            const labelList& fAdd = fAddr();
            return fAdd[fi];
        }

        //- Return start label of this patch in the polyMesh face list
        label zoneId() const
        {
            return zoneId_;
        }

        label& zoneId()
        {
            return zoneId_;
        }

        word indirectPolyPatchType() const
        {
            return indirectPolyPatchType_;
        }

        word& indirectPolyPatchType()
        {
            return indirectPolyPatchType_;
        }

        labelList fAddr() const
        {
            return fAddr_;
        }

        const labelList& fAddr()
        {
            return fAddr_;
        }

        boolList fm() const
        {
            return fm_;
        }

        const boolList& fm()
        {
            return fm_;
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

        //- remove the contributions from the matrix
        virtual void GIBAgglomerate
        (
            scalarField&
        ) const;

        //- Slice list to patch
        template<class T>
        const typename List<T>::subList patchSlice(const UList<T>& l) const
        {
            WarningInFunction
                << "for indirectPolyPatch not supported."
                << endl;
            return typename List<T>::subList(l, this->size(), start_);
        }

        //- Slice Field to patch
        template<class T>
        const typename Field<T>::subField patchSlice(const Field<T>& l) const
        {
            WarningInFunction
                << "for indirectPolyPatch not supported."
                << endl;
            return typename Field<T>::subField(l, this->size(), start_);
        }

        virtual const List<scalar>::subList patchSlice
        (
            const Field<scalar>& l
        ) const
        {
            WarningInFunction
                << "for indirectPolyPatch not supported."
                << endl;
            scalarField fc = scalarField(size());
            forAll(fc, facei)
            {
                label gfI = localToGlobal(facei);
                fc[facei] = l[gfI];
            }
            scalarField::subField sfc = scalarField::subField(fc, this->size(), 0);

            return sfc;
        }
        virtual const List<vector>::subList patchSlice
        (
            const Field<vector>& l
        ) const
        {
            WarningInFunction
                << "for indirectPolyPatch not supported."
                << endl;
            vectorField fc = vectorField(this->size());
            forAll(fc, facei)
            {
                label gfI = localToGlobal(facei);
                fc[facei] = l[gfI];
            }
            vectorField::subField sfc = vectorField::subField(fc, this->size(), 0);

            return sfc;
        }
        virtual const List<scalar>::subList patchSlice
        (
            const List<scalar>& l
        ) const
        {
            return List<scalar>::subList(l, this->size(), 0);
        }
        virtual const List<vector>::subList patchSlice
        (
            const List<vector>& l
        ) const
        {
            return List<vector>::subList(l, this->size(), 0);
        }
        virtual const List<vector> indirectPatchSlice
        (
            const Field<vector>& l
        ) const
        {
            List<vector> fc = List<vector>(this->size());
            forAll(fc, facei)
            {
                label gfI = localToGlobal(facei);
                fc[facei] = l[gfI];
            }

            return fc;
        }
        virtual const List<scalar> indirectPatchSlice
        (
            const Field<scalar>& l
        ) const
        {
            List<scalar> fc = List<scalar>(this->size());
            forAll(fc, facei)
            {
                label gfI = localToGlobal(facei);
                fc[facei] = l[gfI];
            }

            return fc;
        }

        //- Write the indirectPolyPatch data as a dictionary
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

        // Pure functions from indirectPrimitivePatch
            virtual label size() const
            {
                return this->indirectPrimitivePatch::size();
            }
            virtual const labelList& boundaryPoints() const
            {
                return this->indirectPrimitivePatch::boundaryPoints();
            }
            virtual const labelListList& faceFaces() const
            {
                return this->indirectPrimitivePatch::faceFaces();
            }

            virtual const labelListList& edgeFaces() const
            {
                return this->indirectPrimitivePatch::edgeFaces();
            }
            virtual const labelListList& faceEdges() const
            {
                return this->indirectPrimitivePatch::faceEdges();
            }
            virtual const labelListList& pointEdges() const
            {
                return this->indirectPrimitivePatch::pointEdges();
            }
            virtual const labelListList& pointFaces() const
            {
                return this->indirectPrimitivePatch::pointFaces();
            }
            virtual label nPoints() const
            {
                return this->indirectPrimitivePatch::nPoints();
            }
            virtual label nEdges() const
            {
                return this->indirectPrimitivePatch::nEdges();
            }
            virtual const edgeList& edges() const
            {
                return this->indirectPrimitivePatch::edges();
            }
            virtual label nInternalEdges() const
            {
                return this->indirectPrimitivePatch::nInternalEdges();
            }
            virtual bool isInternalEdge(const label edgeI) const
            {
                return this->indirectPrimitivePatch::isInternalEdge(edgeI);
            }
            virtual const List<face>& localFaces() const
            {
                return this->indirectPrimitivePatch::localFaces();
            }
            virtual const labelList& meshPoints() const
            {
                return this->indirectPrimitivePatch::meshPoints();
            }
            virtual const Map<label>& meshPointMap() const
            {
                return this->indirectPrimitivePatch::meshPointMap();
            }
            virtual const Field<point>& localPoints() const
            {
                return this->indirectPrimitivePatch::localPoints();
            }
            const Field<point>& faceNormals() const
            {
                return this->indirectPrimitivePatch::faceNormals();
            }
            const Field<point>& pointNormals() const
            {
                return this->indirectPrimitivePatch::pointNormals();
            }
            const Field<point>& points() const
            {
                return this->indirectPrimitivePatch::points();
            }
            virtual bool empty() const
            {
                return this->indirectPrimitivePatch::empty();
            }
            virtual const labelListList& edgeLoops() const
            {
                return this->indirectPrimitivePatch::edgeLoops();
            }
            virtual const Field<scalar>& faceCurvature() const
            {
                return this->indirectPrimitivePatch::faceCurvature();
            }
            virtual const Field<scalar>& pointCurvature() const
            {
                return this->indirectPrimitivePatch::pointCurvature(fm());
            }

        // Other patch operations

            //- Return label of face in patch from global face label
            //- this is wrong
            virtual label whichFace(const label l) const
            {
                WarningInFunction
                    << "for indirectPolyPatch not supported."
                     << endl;
                return l - start_;
            }

            virtual label whichPoint(const label gp) const
            {
                return this->indirectPrimitivePatch::whichPoint(gp);
            }

            virtual label whichEdge(const edge& ed) const
            {
                return this->indirectPrimitivePatch::whichEdge(ed);
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
        void operator=(const polyPatch&);
        virtual const face& operator[](const label) const;
        virtual face& operator[](const label);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const indirectPolyPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
