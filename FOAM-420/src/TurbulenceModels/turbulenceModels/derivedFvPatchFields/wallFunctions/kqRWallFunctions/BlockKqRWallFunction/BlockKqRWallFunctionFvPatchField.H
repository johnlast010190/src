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

Class
    Foam::BlockKqRWallFunctionFvPatchField

Group
    grpWallFunctions

Description
    This boundary condition provides a suitable condition for turbulence
    \c k, \c q, and \c R fields for the case of high Reynolds number flow using
    wall functions.

    Same as kqRWallFunction. Backward compatibility

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            BlockKqRWallFunction;
    }
    \endverbatim

SeeAlso

SourceFiles
    BlockKqRWallFunctionFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef BlockKqRWallFunctionFvPatchField_H
#define BlockKqRWallFunctionFvPatchField_H

#include "derivedFvPatchFields/wallFunctions/kqRWallFunctions/kqRWallFunction/kqRWallFunctionFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class BlockKqRWallFunctionFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockKqRWallFunctionFvPatchField
:
    public kqRWallFunctionFvPatchField<Type>
{

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("BlockKqRWallFunction");


    // Constructors

        //- Construct from patch and internal field
        BlockKqRWallFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        BlockKqRWallFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // BlockKqRWallFunctionFvPatchField
        //  onto a new patch
        BlockKqRWallFunctionFvPatchField
        (
            const BlockKqRWallFunctionFvPatchField&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        BlockKqRWallFunctionFvPatchField
        (
            const BlockKqRWallFunctionFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new BlockKqRWallFunctionFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        BlockKqRWallFunctionFvPatchField
        (
            const BlockKqRWallFunctionFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new BlockKqRWallFunctionFvPatchField(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
