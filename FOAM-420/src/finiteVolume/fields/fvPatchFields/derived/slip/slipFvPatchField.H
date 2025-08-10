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
    Foam::slipFvPatchField

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition provides a slip constraint.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            slip;
    }
    \endverbatim

SourceFiles
    slipFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef slipFvPatchField_H
#define slipFvPatchField_H

#include "fields/fvPatchFields/basic/basicSymmetry/basicSymmetryFvPatchField.H"
#include "referenceFrames/coordinateFrame.H"
#include "referenceFrames/referenceFrameFvPatch/referenceFrameFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class slipFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class slipFvPatchField
:
    public basicSymmetryFvPatchField<Type>,
    public referenceFrameFvPatch<Type>
{
using referenceFrameFvPatch<Type>::getFrameVelocity;

    //- Private member data

        // Relax
        const scalar relax_;

public:

    //- Runtime type information
    TypeName("slip");


    // Constructors

        //- Construct from patch and internal field
        slipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        slipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given slipFvPatchField onto a new patch
        slipFvPatchField
        (
            const slipFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        slipFvPatchField
        (
            const slipFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const override
        {
            return tmp<fvPatchField<Type>>
            (
                new slipFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        slipFvPatchField
        (
            const slipFvPatchField<Type>&,
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
                new slipFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&) override;

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            ) override;

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&) override;


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const override;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
            ) override;


            //- Return the matrix diagonal coefficients
            virtual tmp<Field<Type>> valueDivInternalCoeffs
            (
                const tmp<scalarField>&
            ) const override;

            //- Return the matrix source coefficients
            virtual tmp<Field<Type>> valueDivBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const override;

        // Attributes

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const override
            {
                return false;
            }

        //- Relax the boundary condition
        virtual void boundaryRelaxMatrix(fvBlockMatrix<Type>&) const override {}

        //- Write
        virtual void write(Ostream&) const override;

    // Member operators

        virtual void operator=(const UList<Type>&) override {}

        virtual void operator=(const fvPatchField<Type>&) override {}
        virtual void operator+=(const fvPatchField<Type>&) override {}
        virtual void operator-=(const fvPatchField<Type>&) override {}
        virtual void operator*=(const fvPatchField<scalar>&) override {}
        virtual void operator/=(const fvPatchField<scalar>&) override {}

        virtual void operator+=(const Field<Type>&) override {}
        virtual void operator-=(const Field<Type>&) override {}

        virtual void operator*=(const Field<scalar>&) override {}
        virtual void operator/=(const Field<scalar>&) override {}

        virtual void operator=(const Type&) override {}
        virtual void operator+=(const Type&) override {}
        virtual void operator-=(const Type&) override {}
        virtual void operator*=(const scalar) override {}
        virtual void operator/=(const scalar) override {}
};

// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<scalarField > slipFvPatchField<scalar>::snGrad() const;

template<>
void slipFvPatchField<scalar>::evaluate
(
    const Pstream::commsTypes commsType
);


template<>
void slipFvPatchField<vector>::boundaryRelaxMatrix
(
    fvBlockMatrix<vector>&
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
