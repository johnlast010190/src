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
    (c) 2020 Esi Ltd.

Class
    Foam::wedgeFvPatchField

Group
    grpConstraintBoundaryConditions

Description
    This boundary condition is similar to the cyclic condition, except that
    it is applied to 2-D geometries.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            wedge;
    }
    \endverbatim

See also
    Foam::cyclicFvPatchField

SourceFiles
    wedgeFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeFvPatchField_H
#define wedgeFvPatchField_H

#include "fields/fvPatchFields/basic/transform/transformFvPatchField.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class wedgeFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class wedgeFvPatchField
:
    public transformFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName(wedgeFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        wedgeFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wedgeFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wedgeFvPatchField onto a new patch
        wedgeFvPatchField
        (
            const wedgeFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wedgeFvPatchField
        (
            const wedgeFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const override
        {
            return tmp<fvPatchField<Type>>
            (
                new wedgeFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wedgeFvPatchField
        (
            const wedgeFvPatchField<Type>&,
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
                new wedgeFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        //- Return gradient at boundary
        virtual tmp<Field<Type>> snGrad() const override;

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        ) override;

        //- Return face-gradient transform diagonal
        virtual tmp<Field<Type>> snGradTransformDiag() const override;

        //- Return patch value for explicit gradient evaluation
        tmp<Field<Type>> gradientBoundaryValue() const override;
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<scalarField > wedgeFvPatchField<scalar>::snGrad() const;

template<>
void wedgeFvPatchField<scalar>::evaluate
(
    const Pstream::commsTypes commsType
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
