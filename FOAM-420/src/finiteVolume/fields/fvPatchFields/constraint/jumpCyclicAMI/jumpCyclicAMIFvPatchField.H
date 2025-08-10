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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::jumpCyclicAMIFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a base class that enforces a cyclic
    condition with a specified 'jump' (or offset) between a pair of boundaries,
    whereby communication between the patches is performed using an arbitrary
    mesh interface (AMI) interpolation.

See also
    Foam::cyclicAMIFvPatchField

SourceFiles
    jumpCyclicAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef jumpCyclicAMIFvPatchField_H
#define jumpCyclicAMIFvPatchField_H

#include "fields/fvPatchFields/constraint/cyclicAMI/cyclicAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class jumpCyclicAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class jumpCyclicAMIFvPatchField
:
    public cyclicAMIFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("jumpCyclicAMI");


    // Constructors

        //- Construct from patch and internal field
        jumpCyclicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        jumpCyclicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given jumpCyclicAMIFvPatchField onto a
        //  new patch
        jumpCyclicAMIFvPatchField
        (
            const jumpCyclicAMIFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        jumpCyclicAMIFvPatchField
        (
            const jumpCyclicAMIFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        jumpCyclicAMIFvPatchField
        (
            const jumpCyclicAMIFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );


    // Member functions

        // Access

            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                return cyclicAMIFvPatchField<Type>::type();
            }

            //- Return local reference cast into the cyclic AMI patch
            virtual const cyclicAMIFvPatch& cyclicAMIPatch() const
            {
                return cyclicAMIFvPatchField<Type>::cyclicAMIPatch();
            }

            //- Return the "jump" across the patch as a "half" field
            virtual tmp<Field<Type>> jump() const = 0;

            //- Return the "jump" across the patch as a "half" field
            virtual Field<Type>& jumpRef() = 0;


        // Evaluation functions

            //- Return neighbour coupled given internal cell data
            virtual tmp<Field<Type>> patchNeighbourField
            (
                const Pstream::commsTypes
            ) const;

            virtual void initInterfaceMatrixUpdate
            (
                const Field<Type>& f1,
                Field<Type>& f2,
                const BlockLduMatrix<Type>& matrix,
                const CoeffField<Type>& coeff,
                const Pstream::commsTypes com
            ) const
            {
                NotImplemented;
            }

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const bool add,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>&,
                const bool add,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const;

};


//- Update result field based on interface functionality
template<>
void jumpCyclicAMIFvPatchField<scalar>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
