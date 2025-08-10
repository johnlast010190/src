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
    Foam::nonConformalErrorFvPatchField

Description
    This boundary condition enforces a non-conformal error condition.

Usage
    Example of the boundary condition specification:

    \verbatim
    <patchName>
    {
        type            nonConformalError;
    }
    \endverbatim

SourceFiles
    nonConformalErrorFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalErrorFvPatchField_H
#define nonConformalErrorFvPatchField_H

#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"
#include "fvMesh/fvPatches/constraint/nonConformalError/nonConformalErrorFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class nonConformalErrorFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalErrorFvPatchField
:
    public zeroGradientFvPatchField<Type>
{
    // Private data

        //- Reference to the fvPatch
        const nonConformalErrorFvPatch& nonConformalErrorFvPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalErrorFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        nonConformalErrorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalErrorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonConformalErrorFvPatchField onto a
        //  new patch
        nonConformalErrorFvPatchField
        (
            const nonConformalErrorFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nonConformalErrorFvPatchField
        (
            const nonConformalErrorFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        nonConformalErrorFvPatchField
        (
            const nonConformalErrorFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchField<Type>>
            (
                new nonConformalErrorFvPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/fvPatchFields/constraint/nonConformalError/nonConformalErrorFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
