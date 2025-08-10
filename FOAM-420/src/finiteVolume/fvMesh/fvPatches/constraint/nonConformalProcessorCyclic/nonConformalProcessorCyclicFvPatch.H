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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::nonConformalProcessorCyclicFvPatch

Description
    Non-conformal processor cyclic FV patch. As nonConformalCyclicFvPatch,
    but the neighbouring patch is on a different processor.

See also
    Foam::nonConformalProcessorCyclicPolyPatch
    Foam::nonConformalCyclicFvPatch

SourceFiles
    nonConformalProcessorCyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicFvPatch_H
#define nonConformalProcessorCyclicFvPatch_H

#include "fvMesh/fvPatches/constraint/processorCyclic/processorCyclicFvPatch.H"
#include "fvMesh/fvPatches/constraint/nonConformalCoupled/nonConformalCoupledFvPatch.H"
#include "nonConformal/polyPatches/nonConformalProcessorCyclic/nonConformalProcessorCyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class nonConformalProcessorCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalProcessorCyclicFvPatch
:
    public processorCyclicFvPatch,
    public nonConformalCoupledFvPatch
{
    // Private data

        //- Reference to the polyPatch
        const nonConformalProcessorCyclicPolyPatch&
            nonConformalProcessorCyclicPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField& w) const override;


public:

    //- Runtime type information
    TypeName(nonConformalProcessorCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        nonConformalProcessorCyclicFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        );


    //- Destructor
    virtual ~nonConformalProcessorCyclicFvPatch();


    // Member Functions

        // Access

            using processorCyclicFvPatch::owner;
            using processorCyclicFvPatch::neighbour;
            using processorCyclicFvPatch::transform;

            //- Poly patch
            const nonConformalProcessorCyclicPolyPatch&
                nonConformalProcessorCyclicPatch() const
            {
                return nonConformalProcessorCyclicPolyPatch_;
            }

            //- Is this patch coupled? For NCC patches, the poly mesh is
            //  considered non-coupled whilst the finite volume mesh is
            //  considered coupled.
            virtual bool coupled() const override
            {
                return Pstream::parRun();
            }

            //- Return the start label of this patch in the polyMesh face list.
            //  Raise an error unless the patch is empty.
            virtual label start() const override
            {
                return nonConformalFvPatch::start();
            }

            //- Return the size
            virtual label size() const override
            {
                return nonConformalFvPatch::size();
            }

            //- Return the face-cells
            virtual const labelUList& faceCells() const override
            {
                return nonConformalFvPatch::faceCells();
            }


        //- Return delta (P to N) vectors across coupled patch
        virtual tmp<vectorField> delta() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
