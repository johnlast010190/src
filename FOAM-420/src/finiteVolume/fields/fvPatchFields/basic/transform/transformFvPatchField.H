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
    (c) 2019 Esi Ltd.

Class
    Foam::transformFvPatchField

Group
    grpGenericBoundaryConditions

Description
    Foam::transformFvPatchField

SourceFiles
    transformFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef transformFvPatchField_H
#define transformFvPatchField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class transformFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class transformFvPatchField
:
    public fvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("transform");


    // Constructors

        //- Construct from patch and internal field
        transformFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        transformFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given transformFvPatchField<Type>
        //  onto a new patch
        transformFvPatchField
        (
            const transformFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        transformFvPatchField
        (
            const transformFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const = 0;

        //- Construct as copy setting internal field reference
        transformFvPatchField
        (
            const transformFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>&
        ) const = 0;


    // Member functions

        // Attributes

            //- Return true if the value of the patch field
            //  is altered by assignment
            virtual bool assignable() const
            {
                return true;
            }


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const = 0;

            //- Return face-gradient transform diagonal
            virtual tmp<Field<Type>> snGradTransformDiag() const = 0;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;

            virtual tmp<CoeffField<Type>> gradientInternalBCoeffs() const
            {
                return this->fvPatchField<Type>::gradientInternalBCoeffs();
            };

            virtual tmp<Field<Type>> gradientBoundaryBCoeffs() const
            {
                return this->fvPatchField<Type>::gradientBoundaryBCoeffs();
            };


    // Member operators

        virtual void operator=(const fvPatchField<Type>&);
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<scalarField> transformFvPatchField<scalar>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const;

template<>
tmp<scalarField> transformFvPatchField<scalar>::gradientInternalCoeffs() const;


template<>
tmp<CoeffField<vector>> transformFvPatchField<vector>::gradientInternalBCoeffs
(
) const;


template<>
tmp<Field<vector>> transformFvPatchField<vector>::gradientBoundaryBCoeffs
(
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
