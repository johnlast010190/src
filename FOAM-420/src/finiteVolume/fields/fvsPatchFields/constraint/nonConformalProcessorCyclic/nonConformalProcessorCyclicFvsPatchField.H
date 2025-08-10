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
    Foam::nonConformalProcessorCyclicFvsPatchField

Description
    Foam::nonConformalProcessorCyclicFvsPatchField

SourceFiles
    nonConformalProcessorCyclicFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicFvsPatchField_H
#define nonConformalProcessorCyclicFvsPatchField_H

#include "fields/fvsPatchFields/constraint/processorCyclic/processorCyclicFvsPatchField.H"
#include "fvMesh/fvPatches/constraint/nonConformalProcessorCyclic/nonConformalProcessorCyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class nonConformalProcessorCyclicFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalProcessorCyclicFvsPatchField
:
    public processorCyclicFvsPatchField<Type>
{
    // Private data

        //- Local reference cast into the processor patch
        const nonConformalProcessorCyclicFvPatch& procPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalProcessorCyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        nonConformalProcessorCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch and internal field and patch field
        nonConformalProcessorCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalProcessorCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonConformalProcessorCyclicFvsPatchField
        //  onto a new patch
        nonConformalProcessorCyclicFvsPatchField
        (
            const nonConformalProcessorCyclicFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nonConformalProcessorCyclicFvsPatchField
        (
            const nonConformalProcessorCyclicFvsPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        nonConformalProcessorCyclicFvsPatchField
        (
            const nonConformalProcessorCyclicFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const override
        {
            return tmp<fvsPatchField<Type>>
            (
                new nonConformalProcessorCyclicFvsPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~nonConformalProcessorCyclicFvsPatchField();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/fvsPatchFields/constraint/nonConformalProcessorCyclic/nonConformalProcessorCyclicFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
