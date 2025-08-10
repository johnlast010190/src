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
    (c) 2015 Esi Ltd.

Class
    Foam::outletFixedValueFvPatchField

Group
    grpOutletBoundaryConditions

Description
    This boundary condition provides a weak fixed value outflow condition,
    with strong boundary conditions being applied in the case of inflow

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        value        | Patch face values       | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            outletFixedValue;
        phi             phi;
        value           uniform (1 0 0);
    }
    \endverbatim

    The mode of operation (strong/weak) is determined by the sign of the
    flux across the patch faces.

Note
    Sign conventions:
    - positive flux (out of domain): apply zero-gradient on convection condition
    - negative flux (into of domain): apply the fixed value on convection

SourceFiles
    outletFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef outletFixedValueFvPatchField_H
#define outletFixedValueFvPatchField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class outletFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class outletFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        word phiName_;


public:

    //- Runtime type information
    TypeName("outletFixedValue");


    // Constructors

        //- Construct from patch and internal field
        outletFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outletFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given outletFixedValueFvPatchField
        //  onto a new patch
        outletFixedValueFvPatchField
        (
            const outletFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outletFixedValueFvPatchField
        (
            const outletFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new outletFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        outletFixedValueFvPatchField
        (
            const outletFixedValueFvPatchField<Type>&,
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
                new outletFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions


        // Evaluation functions


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



        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
