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
    (c) 2017 Esi Ltd.

Class
    Foam::surfaceNormalFixedValueFvPatchVectorField

Group
    grpGenericBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a surface-normal vector boundary condition
    by its magnitude. Optionally, a time-varying uniform magnitude may be
    specified via a Function1.

Usage
    \table
        Property     | Description             | Required    | Default value
        refValue     | reference value         | yes         |
        profile      | uniform function        | no          | null
    \endtable

    Example of the constant boundary condition specification:
    \verbatim
    <patchName>
    {
        type            surfaceNormalFixedValue;
        refValue        uniform -10;           // 10 INTO the domain
    }
    \endverbatim


    Example of the time-varying boundary condition specification:
    \verbatim
    <patchName>
    {
        type            surfaceNormalFixedValue;
        refValue        uniform -10;           // 10 INTO the domain
        profile         table
        (
            (0 -10)
            (1 -20)
            (2 -40)
        );
    }
    \endverbatim

Note
    Sign conventions:
    - the value is positive for outward-pointing vectors

See also
    Foam::fixedValueFvPatchField

SourceFiles
    surfaceNormalFixedValueFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceNormalFixedValueFvPatchVectorField_H
#define surfaceNormalFixedValueFvPatchVectorField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "primitives/bools/Switch/Switch.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class surfaceNormalFixedValueFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class surfaceNormalFixedValueFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        scalarField refValue_;

        //- time varying Function1, optional
        autoPtr<Function1<scalar>> profile_;

        // enable weak formulation for outlets
        Switch outlet_;


public:

    //- Runtime type information
    TypeName("surfaceNormalFixedValue");


    // Constructors

        //- Construct from patch and internal field
        surfaceNormalFixedValueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceNormalFixedValueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  surfaceNormalFixedValueFvPatchVectorField
        //  onto a new patch
        surfaceNormalFixedValueFvPatchVectorField
        (
            const surfaceNormalFixedValueFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        surfaceNormalFixedValueFvPatchVectorField
        (
            const surfaceNormalFixedValueFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const override
        {
            return tmp<fvPatchVectorField>
            (
                new surfaceNormalFixedValueFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        surfaceNormalFixedValueFvPatchVectorField
        (
            const surfaceNormalFixedValueFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchVectorField>
            (
                new surfaceNormalFixedValueFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            ) override;

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            ) override;

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&) override;


        // Evaluation functions

            //- Is defined in frame supported?
            virtual bool definedInFrameSupported() const override
            {
                return true;
            }

            //- The velocity normal to the patch should no longer be rotated
            // the boundary itself is already localised
            virtual bool isDefinedInFrame() const override
            {
                return false;
            }

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs() override;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<vectorField> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const override;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<vectorField> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const override;

            //- Return the matrix diagonal coefficients
            virtual tmp<Field<vector>> valueDivInternalCoeffs
            (
                const tmp<scalarField>&
            ) const override;

            //- Return the matrix source coefficients
            virtual tmp<Field<vector>> valueDivBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const override;

        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
