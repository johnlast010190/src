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
    (c) 2018 Esi Ltd.
    (c) 2016-2017 Wikki Ltd

Class
    Foam::coupledFaPatch

Author
    Zeljko Tukovic and Hrvoje Jasak

Description
    coupledFaPatch is an abstract base class for patches that couple regions
    of the computational domain e.g. cyclic, arbitrary interfaces, sliding
    interfaces and processor-processor links.

SourceFiles
    coupledFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFaPatch_H
#define coupledFaPatch_H

#include "matrices/lduMatrix/lduAddressing/lduInterface/lduInterface.H"
#include "faMesh/faPatches/faPatch/faPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coupledFaPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledFaPatch
:
    public lduInterface,
    public faPatch
{
    // Private data

        //- offset (distance) vector from one side of the couple to the other
        mutable vectorField separation_;

        //- Face transformation tensor
        mutable tensorField forwardT_;

        //- Neighbour-cell transformation tensor
        mutable tensorField reverseT_;


protected:

    // Protected Member Functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const = 0;

        //- Calculate the uniform transformation tensors
        void calcTransformTensors
        (
            const vector& Cf,
            const vector& Cr,
            const vector& nf,
            const vector& nr
        ) const;

        //- Calculate the transformation tensors
        void calcTransformTensors
        (
            const vectorField& Cf,
            const vectorField& Cr,
            const vectorField& nf,
            const vectorField& nr
        ) const;


public:

    //- Runtime type information
    TypeName("coupled");


    // Constructors

        //- Construct from components
        coupledFaPatch
        (
            const word& name,
            const labelList& edgeLabels,
            const label index,
            const faBoundaryMesh& bm,
            const label ngbPolyPatchIndex
        )
        :
            faPatch(name, edgeLabels, index, bm, ngbPolyPatchIndex)
        {}

        //- Construct from dictionary
        coupledFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        )
        :
            faPatch(name, dict, index, bm)
        {}


    // Destructor

        virtual ~coupledFaPatch();


    // Member Functions

        // Access

            //- Return true because this patch is coupled
            virtual bool coupled() const
            {
                return true;
            }


            //- Are the coupled planes separated
            bool separated() const
            {
                return separation_.size();
            }

            //- Return the offset (distance) vector from one side of the couple
            //  to the other
            const vectorField& separation() const
            {
                if (!separation_.size())
                {
                    FatalErrorInFunction
                        << "Coupled patches are not separated"
                        << abort(FatalError);
                }

                return separation_;
            }

            //- Return face transformation tensor
            const tensorField& forwardT() const
            {
                if (!forwardT_.size())
                {
                    FatalErrorInFunction
                        << "Coupled planes do not need transformation"
                        << abort(FatalError);
                }

                return forwardT_;
            }

            //- Return neighbour-cell transformation tensor
            const tensorField& reverseT() const
            {
                if (!reverseT_.size())
                {
                    FatalErrorInFunction
                        << "Coupled planes do not need transformation"
                        << abort(FatalError);
                }

                return reverseT_;
            }

            //- Are the cyclic planes parallel
            bool parallel() const
            {
                return forwardT_.size() == 0;
            }


        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&) = 0;

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&) = 0;

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&) = 0;

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&) = 0;


        // Access functions for demand driven data

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const = 0;


        // Interface transfer functions

            //- Return faceCell addressing: lduInterface virtual function
            virtual const labelUList& faceCells() const
            {
                return edgeFaces();
            }

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const = 0;

            //- Inherit initInternalFieldTransfer from lduInterface
             using lduInterface::initInternalFieldTransfer;

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& interfaceData
            ) const = 0;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                labelUList& iF
            ) const
            {}

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
