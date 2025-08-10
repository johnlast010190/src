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
    (c) 2016-2017 Wikki Ltd
    (c) 2018-2022 Esi Ltd.

Class
    Foam::processorFaPatch

Description
    Processor patch.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    processorFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorFaPatch_H
#define processorFaPatch_H

#include "faMesh/faPatches/basic/coupled/coupledFaPatch.H"
#include "matrices/lduMatrix/lduAddressing/lduInterface/processorLduInterface.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class processorFaPatch Declaration
\*---------------------------------------------------------------------------*/

class processorFaPatch
:
    public coupledFaPatch,
    public processorLduInterface
{
    // Private data

        int myProcNo_;
        int neighbProcNo_;

        //- Processor-neighbour patch edge centres
        vectorField nbrEdgeCentres_;

        //- Processor-neighbour patch edge lengths
        vectorField nbrEdgeLengths_;

        //- Processor-neighbour patch neighbour face centres
        vectorField nbrEdgeFaceCentres_;

        //- Corresponding neighbouring local point label for every local point
        //  (so localPoints()[i] == neighb.localPoints()[nbrPoints_[i]])
        mutable labelList* nbrPointsPtr_;

        mutable labelList* nbrEdgesPtr_;

        //- The set of labels of the processor patch points which are
        //  non-global, i.e. present in this processor patch
        mutable labelList* nonGlobalPatchPointsPtr_;

protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

        //- Make patch face - neighbour cell distances
        void makeDeltaCoeffs(scalarField&) const;

        //- Find non-globa patch points
        void makeNonGlobalPatchPoints() const;

        //- Initialise the calculation of the patch geometry
        void initGeometry(PstreamBuffers&);

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
//     TypeName(processorPolyPatch::typeName_());
    TypeName("processor");


    // Constructors

        //- Construct from components
        processorFaPatch
        (
            const word& name,
            const labelList& edgeLabels,
            const label index,
            const faBoundaryMesh& bm,
            const label ngbPolyPatchIndex,
            const label myProcNo,
            const label neighbProcNo
        )
        :
            coupledFaPatch(name, edgeLabels, index, bm, ngbPolyPatchIndex),
            myProcNo_(myProcNo),
            neighbProcNo_(neighbProcNo),
            nbrEdgeCentres_(),
            nbrEdgeLengths_(),
            nbrEdgeFaceCentres_(),
            nbrPointsPtr_(nullptr),
            nbrEdgesPtr_(nullptr),
            nonGlobalPatchPointsPtr_(nullptr)
        {}

        //- Construct from dictionary
        processorFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        )
        :
            coupledFaPatch(name, dict, index, bm),
            myProcNo_(readLabel(dict.lookup("myProcNo"))),
            neighbProcNo_(readLabel(dict.lookup("neighbProcNo"))),
            nbrEdgeCentres_(),
            nbrEdgeLengths_(),
            nbrEdgeFaceCentres_(),
            nbrPointsPtr_(nullptr),
            nbrEdgesPtr_(nullptr),
            nonGlobalPatchPointsPtr_(nullptr)
        {}

    // Destructor

        virtual ~processorFaPatch();


    // Member functions

        virtual label comm() const
        {
            return UPstream::worldComm;
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

        //- Return message tag used for sending
        virtual int tag() const
        {
            // Use a different tag to protect communication that is split
            // between initEvaluate and evaluate
            return UPstream::msgType()+2;
        }

        //- Return true if running parallel
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

        //- Return face transformation tensor
        virtual const tensorField& forwardT() const
        {
            return coupledFaPatch::forwardT();
        }

        //- Return dummy null transform
        virtual const transformer& transform() const
        {
            return transformer::null;
        }

        //- Return delta (P to N) vectors across coupled patch
        virtual tmp<vectorField> delta() const;

        //- Return processor-neighbour patch edge centres
        const vectorField& nbrEdgeCentres() const
        {
            return nbrEdgeCentres_;
        }

        //- Return processor-neighbour patch edge lengths
        const vectorField& nbrEdgeLengths() const
        {
            return nbrEdgeLengths_;
        }

        //- Return processor-neighbour patch neighbour face centres
        const vectorField& nbrEdgeFaceCentres() const
        {
            return nbrEdgeFaceCentres_;
        }

        //- Return neighbour point labels. This is for my local point the
        //  corresponding local point on the other side. Call
        //  faBoundaryMesh::updateMesh() on all processors
        //  before using this.
        const labelList& nbrPoints() const;

        const labelList& nbrEdges() const;

        //- Return the set of labels of the processor patch points which are
        //  non-global, i.e. present in this processorFaPatch
        const labelList& nonGlobalPatchPoints() const;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Initialise interface data transfer
            virtual void initTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& interfaceData
            ) const;

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& interfaceData
            ) const;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;

            //- Write the patch data as a dictionary
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
