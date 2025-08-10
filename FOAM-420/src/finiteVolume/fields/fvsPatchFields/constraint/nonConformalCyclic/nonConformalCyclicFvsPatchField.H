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
    Foam::nonConformalCyclicFvsPatchField

Description
    Constraint fvsPatchField for nonConformalCyclic couplings.

SourceFiles
    nonConformalCyclicFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCyclicFvsPatchField_H
#define nonConformalCyclicFvsPatchField_H

#include "fields/fvsPatchFields/constraint/cyclic/cyclicFvsPatchField.H"
#include "fvMesh/fvPatches/constraint/nonConformalCyclic/nonConformalCyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class nonConformalCyclicFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalCyclicFvsPatchField
:
    public cyclicFvsPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(nonConformalCyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and nonConformalCyclic field
        nonConformalCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, nonConformalCyclic field and dictionary
        nonConformalCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonConformalCyclicFvsPatchField onto a
        //  new patch
        nonConformalCyclicFvsPatchField
        (
            const nonConformalCyclicFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting nonConformalCyclic field reference
        nonConformalCyclicFvsPatchField
        (
            const nonConformalCyclicFvsPatchField<Type>&
        ) = delete;

        //- Copy constructor setting nonConformalCyclic field reference
        nonConformalCyclicFvsPatchField
        (
            const nonConformalCyclicFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting nonConformalCyclic field
        //  reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const override
        {
            return tmp<fvsPatchField<Type>>
            (
                new nonConformalCyclicFvsPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/fvsPatchFields/constraint/nonConformalCyclic/nonConformalCyclicFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
