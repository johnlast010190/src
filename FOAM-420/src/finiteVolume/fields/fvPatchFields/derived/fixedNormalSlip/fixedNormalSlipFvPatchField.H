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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::fixedNormalSlipFvPatchField

Group
    grpGenericBoundaryConditions grpWallBoundaryConditions

Description
    This boundary condition sets the patch-normal component to a fixed value.

Usage
    \table
        Property     | Description             | Required    | Default value
        fixedValue   | fixed value             | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedNormalSlip;
        fixedValue      uniform 0;     // example entry for a scalar field
    }
    \endverbatim

See also
    Foam::transformFvPatchField

SourceFiles
    fixedNormalSlipFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedNormalSlipFvPatchField_H
#define fixedNormalSlipFvPatchField_H

#include "fields/fvPatchFields/basic/transform/transformFvPatchField.H"
#include "referenceFrames/coordinateFrame.H"
#include "referenceFrames/referenceFrameFvPatch/referenceFrameFvPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class fixedNormalSlipFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedNormalSlipFvPatchField
:
    public transformFvPatchField<Type>,
    public referenceFrameFvPatch<Type>
{
using referenceFrameFvPatch<Type>::getFrameVelocity;

protected:

        //- Value the normal component of which the boundary is set to
        Field<Type> fixedValue_;

private:

    //- Private member data

        // Relax
        const scalar relax_;

public:

    //- Runtime type information
    TypeName("fixedNormalSlip");


    // Constructors

        //- Construct from patch and internal field
        fixedNormalSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedNormalSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedNormalSlipFvPatchField
        //  onto a new patch
        fixedNormalSlipFvPatchField
        (
            const fixedNormalSlipFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedNormalSlipFvPatchField
        (
            const fixedNormalSlipFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const override
        {
            return tmp<fvPatchField<Type>>
            (
                new fixedNormalSlipFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedNormalSlipFvPatchField
        (
            const fixedNormalSlipFvPatchField<Type>&,
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
                new fixedNormalSlipFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const override
            {
                return false;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            ) override;

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            ) override;

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&) override;


        // Return defining fields

            virtual Field<Type>& fixedValue()
            {
                return fixedValue_;
            }

            virtual const Field<Type>& fixedValue() const
            {
                return fixedValue_;
            }


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const override;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
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

            //- Return face-gradient transform diagonal
            virtual tmp<Field<Type>> snGradTransformDiag() const override;

        //- Relax the boundary condition
        virtual void boundaryRelaxMatrix(fvBlockMatrix<Type>&) const override {}

        //- Write
        virtual void write(Ostream&) const override;

        //- Is defined in frame supported?
        virtual bool definedInFrameSupported() const override
        {
            return true;
        }


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
void fixedNormalSlipFvPatchField<vector>::boundaryRelaxMatrix
(
    fvBlockMatrix<vector>&
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
