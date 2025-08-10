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
    Foam::fixedJumpAMIFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a jump condition, across non-conformal
    cyclic path-pairs, employing an arbitraryMeshInterface (AMI).

    The jump is specified as a fixed value field, applied as an offset to the
    'owner' patch.

Usage
    \table
        Property     | Description             | Required    | Default value
        patchType    | underlying patch type should be \c cyclic| yes |
        jump         | current jump value      | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedJumpAMI;
        patchType       cyclic;
        jump            uniform 10;
    }
    \endverbatim

    The above example shows the use of a fixed jump of '10'.

Note
     The underlying \c patchType should be set to \c cyclicAMI

See also
    Foam::jumpCyclicAMIFvPatchField

SourceFiles
    fixedJumpAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedJumpAMIFvPatchField_H
#define fixedJumpAMIFvPatchField_H

#include "fields/fvPatchFields/constraint/jumpCyclicAMI/jumpCyclicAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fixedJumpAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedJumpAMIFvPatchField
:
    public jumpCyclicAMIFvPatchField<Type>
{

protected:

    // Protected data

        //- "jump" field
        Field<Type> jump_;


public:

    //- Runtime type information
    TypeName("fixedJumpAMI");


    // Constructors

        //- Construct from patch and internal field
        fixedJumpAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedJumpAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedJumpAMIFvPatchField onto a
        //  new patch
        fixedJumpAMIFvPatchField
        (
            const fixedJumpAMIFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedJumpAMIFvPatchField
        (
            const fixedJumpAMIFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new fixedJumpAMIFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedJumpAMIFvPatchField
        (
            const fixedJumpAMIFvPatchField<Type>&,
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
                new fixedJumpAMIFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the "jump" across the patch
            virtual tmp<Field<Type>> jump() const;

            virtual Field<Type>& jumpRef()
            {
                return jump_;
            };

            virtual tmp<Field<Type>> patchNeighbourField
            (
                const Pstream::commsTypes commsType
            ) const
            {
                return jumpCyclicAMIFvPatchField<Type>::patchNeighbourField
                (
                    commsType
                );
            }

            //- Return internal field next to patch as patch field
            virtual tmp<Field<Type>> patchInternalField() const
            {
                return cyclicAMIFvPatchField<Type>::patchInternalField();
            }

            virtual void updateInterfaceMatrix
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
            ) const
            {
                NotImplemented;
            }

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

            virtual void addToInternalCoupledField
            (
                Field<Type>& f1,
                const Field<Type>& f2,
                const CoeffField<Type>& coeff
            ) const
            {
                NotImplemented;
            }

            virtual void addToInternalCoupledCoveredField
            (
                Field<Type>& f1,
                const Field<Type>& f2,
                const CoeffField<Type>& coeff
            ) const
            {
                NotImplemented;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchField<Type>&, const labelList&);

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
