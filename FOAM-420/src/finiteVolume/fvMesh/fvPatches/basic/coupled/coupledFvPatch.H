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
    (c) 2022 Esi Ltd.

Class
    Foam::coupledFvPatch

Description
    An abstract base class for patches that couple regions of the
    computational domain e.g. cyclic and processor-processor links.

SourceFiles
    coupledFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFvPatch_H
#define coupledFvPatch_H

#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
#include "matrices/lduMatrix/lduAddressing/lduInterface/lduInterface.H"
#include "meshes/polyMesh/polyPatches/basic/coupled/coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coupledFvPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledFvPatch
:
    public lduInterface,
    public fvPatch
{
    // Private data

        const coupledPolyPatch& coupledPolyPatch_;


protected:

    // Protected Member Functions

        //- Make patch weighting factors
        void makeWeights
        (
            scalarField& w,
            const vectorField& nbrSf,
            const vectorField& nbrDelta
        ) const;

        //- Make patch weighting factors
        virtual void makeWeights(scalarField& w) const = 0;

        //- Return delta (P to N) vectors across coupled patch
        tmp<vectorField> delta
        (
            const vectorField& nbrDelta
        ) const;


public:

    friend class surfaceInterpolation;


    //- Runtime type information
    TypeName(coupledPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        coupledFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            coupledPolyPatch_(refCast<const coupledPolyPatch>(patch))
        {}


    //- Destructor
    virtual ~coupledFvPatch();


    // Member Functions

        // Access

            //- Return true because this patch is coupled
            virtual bool coupled() const
            {
                return coupledPolyPatch_.coupled();
            }

            //- Does this side own the patch?
            virtual bool owner() const
            {
                return coupledPolyPatch_.owner();
            }

            //- Does the coupled side own the patch?
            virtual bool neighbour() const
            {
                return !owner();
            }

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const = 0;

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return fvPatch::faceCells();
            }

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const = 0;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const = 0;

            //- Inherit initInternalFieldTransfer from lduInterface
            using lduInterface::initInternalFieldTransfer;

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
