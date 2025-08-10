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
    (c) 2022-2023 Esi Ltd.

Class
    Foam::nonConformalProcessorCyclicFvPatchField

Description
    This boundary condition enables processor communication across cyclic
    patches.

Usage
    Example of the boundary condition specification:

    \verbatim
    <patchName>
    {
        type            nonConformalProcessorCyclic;
    }
    \endverbatim

See also
    Foam::processorCyclicFvPatchField

SourceFiles
    nonConformalProcessorCyclicFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicFvPatchField_H
#define nonConformalProcessorCyclicFvPatchField_H

#include "fields/fvPatchFields/constraint/processorCyclic/processorCyclicFvPatchField.H"
#include "fvMesh/fvPatches/constraint/nonConformalProcessorCyclic/nonConformalProcessorCyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class nonConformalProcessorCyclicFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalProcessorCyclicFvPatchField
:
    public processorCyclicFvPatchField<Type>
{
    // Private data

        //- Local reference cast into the processor patch
        const nonConformalProcessorCyclicFvPatch& procPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalProcessorCyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        nonConformalProcessorCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        nonConformalProcessorCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalProcessorCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonConformalProcessorCyclicFvPatchField
        //  onto a new patch
        nonConformalProcessorCyclicFvPatchField
        (
            const nonConformalProcessorCyclicFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nonConformalProcessorCyclicFvPatchField
        (
            const nonConformalProcessorCyclicFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        nonConformalProcessorCyclicFvPatchField
        (
            const nonConformalProcessorCyclicFvPatchField<Type>&,
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
                new nonConformalProcessorCyclicFvPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~nonConformalProcessorCyclicFvPatchField();


    // Member Functions

        // Access

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const override
            {
                return procPatch_.transform();
            }

            //- Return rank of component for transform
            virtual int rank() const override
            {
                return pTraits<Type>::rank;
            }


        // Block coupled interface functionality

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
            ) const override
            {
                NotImplemented;
            }

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
            ) const override
            {
                NotImplemented;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformalProcessorCyclicFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
