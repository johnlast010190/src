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
    (c) 2015 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::processorPolyPatch

Description
    Neighbour processor patch.

    Note: morph patch face ordering tries to do a geometric ordering.
    (assumes faces coincident) Hence will have problems when cyclics
    are present.

SourceFiles
    processorPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorPolyPatch_H
#define processorPolyPatch_H

#include "meshes/polyMesh/polyPatches/basic/coupled/coupledPolyPatch.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMesh.H"
#include "meshes/meshShapes/face/faceListFwd.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class processorPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class processorPolyPatch
:
    public coupledPolyPatch
{
    // Private data

        int myProcNo_;
        int neighbProcNo_;
        label neighbProPatchID_;

        //- Processor-neighbour patch face centres
        vectorField nbrFaceCentres_;

        //- Processor-neighbour patch face areas
        vectorField nbrFaceAreas_;

        //- Processor-neighbour patch neighbour cell centres
        vectorField nbrFaceCellCentres_;

        //- Corresponding neighbouring local point label for every local point
        //  (so localPoints()[i] == neighb.localPoints()[nbrPoints_[i]])
        mutable autoPtr<labelList> nbrPointsPtr_;

        //- Corresponding neighbouring local edge label for every local edge
        //  (so edges()[i] == neighb.edges()[nbrEdges_[i]])
        mutable autoPtr<labelList> nbrEdgesPtr_;

protected:

    // Protected constructors

        //- Construct from components with specified name
        processorPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const int myProcNo,
            const int neighbProcNo,
            const word& patchType = typeName
        );


    // Protected Member functions

        //- Initialise the calculation of the patch geometry
        void initCalcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("processor");


    // Constructors

        //- Construct from components with automatically generated standard name
        processorPolyPatch
        (
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const int myProcNo,
            const int neighbProcNo,
            const word& patchType = typeName
        );

        //- Construct from dictionary
        processorPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        processorPolyPatch(const processorPolyPatch&, const polyBoundaryMesh&);

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        processorPolyPatch
        (
            const processorPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct given the original patch and a map
        processorPolyPatch
        (
            const processorPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new processorPolyPatch(*this, bm));
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
                new processorPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
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
                new processorPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~processorPolyPatch();


    // Member functions

        //- Return true only if this is a parallel run
        virtual bool coupled() const
        {
            if (Pstream::parRun())
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        //- Return processor number
        int myProcNo() const
        {
            return myProcNo_;
        }

        //- Return neigbour processor number
        int neighbProcNo() const
        {
            return neighbProcNo_;
        }

        //- Return label of the patch in the neighbour processor
        label& neighbProPatchID()
        {
            return neighbProPatchID_;
        }

        //- Does the processor own the patch ?
        virtual bool owner() const
        {
            return (myProcNo_ < neighbProcNo_);
        }

        //- Is the processor the patch neighbour ?
        bool neighbour() const
        {
            return !owner();
        }

        //- Return the name of a processorPolyPatch
        //  constructed from the processor IDs
        static word newName
        (
            const label myProcNo,
            const label neighbProcNo
        );

        //- Return processor-neighbour patch face centres
        const vectorField& nbrFaceCentres() const
        {
            return nbrFaceCentres_;
        }

        //- Return processor-neighbour patch face areas
        const vectorField& nbrFaceAreas() const
        {
            return nbrFaceAreas_;
        }

        //- Return processor-neighbour patch neighbour cell centres
        const vectorField& nbrFaceCellCentres() const
        {
            return nbrFaceCellCentres_;
        }

        //- Return neighbour point labels. WIP.
        const labelList& nbrPoints() const;

        //- Return neighbour edge labels. WIP.
        const labelList& nbrEdges() const;

        //- Return message tag to use for communication
        virtual int tag() const
        {
            return Pstream::msgType();
        }

        //- Return communicator used for communication
        virtual label comm() const
        {
            return boundaryMesh().mesh().comm();
        }

        //- Return null transform between processor patches
        virtual const transformer& transform() const
        {
            return transformer::null;
        }

        //- Get the number of vertices face f needs to be rotated such that
        //  its f[0] point aligns with given anchor (within tol).
        //  Legacy method for old geometric ordering
        static label getRotation
        (
            const pointField& points,
            const face& f,
            const point& anchor,
            const scalar tol
        );

        //- Use old geometric match when detect issue with topological match.
        //  Legacy method for old geometric ordering
        bool geometricMatch
        (
            const primitivePatch& pp,
            const vectorField& masterCtrs,
            const vectorField& masterNormals,
            const vectorField& masterAnchors,
            const vectorField& masterFacePointAverages,
            labelList& faceMap,
            labelList& rotation
        ) const;

        //- Initialise ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder(PstreamBuffers&, const primitivePatch&) const;

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
        ) const;


        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
