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

Class
    Foam::indirectPatchFvPatch

Description
    Foam::indirectPatchFvPatch

SourceFiles
    indirectPatchFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef indirectPatchFvPatch_H
#define indirectPatchFvPatch_H

#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"
#include "matrices/lduMatrix/lduAddressing/lduInterface/lduInterface.H"
#include "fvMesh/fvBoundaryMesh/fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class indirectPatchFvPatch Declaration
\*---------------------------------------------------------------------------*/

class indirectPatchFvPatch
:
    public fvPatch,
    public lduInterface
{

        const indirectPolyPatch& indirectPolyPatch_;

public:

    //- Runtime type information
    TypeName(indirectPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        indirectPatchFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            indirectPolyPatch_(refCast<const indirectPolyPatch>(patch))
        {}

    // Member Functions

        // Access

            label nbrPatchID() const
            {
                return indirectPolyPatch_.nbrPatchID();
            }

            const indirectPatchFvPatch& neighbFvPatch() const
            {
                return refCast<const indirectPatchFvPatch>
                (
                    this->boundaryMesh()[indirectPolyPatch_.nbrPatchID()]
                );
            }

            // functions of lduInterface class
            virtual const labelUList& faceCells() const
            {
                return this->patch().faceCells();
            }

        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const
            {
                return patchInternalField(internalData);
            }

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const
            {
                return neighbFvPatch().patchInternalField(iF);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
