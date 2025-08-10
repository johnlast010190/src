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
    (c) 2011 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::processorCyclicFvPatch

Description
    Processor patch.

SourceFiles
    processorCyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicFvPatch_H
#define processorCyclicFvPatch_H

#include "meshes/polyMesh/polyPatches/constraint/processorCyclic/processorCyclicPolyPatch.H"
#include "fvMesh/fvPatches/constraint/processor/processorFvPatch.H"
#include "fvMesh/fvPatches/constraint/cyclic/cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class processorCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class processorCyclicFvPatch
:
    public processorFvPatch
{
    // Private Data

        const processorCyclicPolyPatch& procPolyPatch_;

public:

    //- Runtime type information
    TypeName(processorCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        processorCyclicFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        )
        :
            processorFvPatch(patch, bm),
            procPolyPatch_(refCast<const processorCyclicPolyPatch>(patch))
        {}


    // Member functions

        //- Return message tag used for sending
        virtual int tag() const
        {
            // Allocate from Pstream?
            return procPolyPatch_.tag();
        }

        //- Return the referring patch ID
        label referPatchID() const
        {
            return procPolyPatch_.referPatchID();
        }

        //- Return a reference to the referring patch
        const cyclicFvPatch& referPatch() const
        {
            const fvPatch& fvp = this->boundaryMesh()[referPatchID()];
            return refCast<const cyclicFvPatch>(fvp);
        }

        //- Return transformation between the coupled patches
        virtual const transformer& transform() const
        {
            return procPolyPatch_.transform();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
