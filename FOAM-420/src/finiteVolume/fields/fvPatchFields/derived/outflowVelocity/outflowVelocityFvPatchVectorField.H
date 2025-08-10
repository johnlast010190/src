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
    (c) 2021 Esi Ltd.

Class
    Foam::outflowVelocityFvPatchVectorField

Description
    Outflow BC for velocity. Should be used with zeroGradient for p.
    It adjust(scales) normal component to satisfy continuity.
    Zero diffussion flux is imposed. Negative components are cropped.

SourceFiles
    outflowVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef outflowVelocityFvPatchVectorField_H
#define outflowVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class outflowVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class outflowVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- patch relax in block solver
        scalar relax_;

        //- scale normal to satify continuity
        scalar scale_;

        //-
        scalar massCorr_;


public:

   //- Runtime type information
   TypeName("outflowVelocity");


   // Constructors

        //- Construct from patch and internal field
        outflowVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outflowVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  outflowVelocityFvPatchVectorField
        //  onto a new patch
        outflowVelocityFvPatchVectorField
        (
            const outflowVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outflowVelocityFvPatchVectorField
        (
            const outflowVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const override
        {
            return tmp<fvPatchVectorField>
            (
                new outflowVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        outflowVelocityFvPatchVectorField
        (
            const outflowVelocityFvPatchVectorField&,
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
                new outflowVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions
        // Attributes

            virtual bool fixesValue() const override
            {
                return false;
            }

        // Evaluation functions

            void scaleProfile(const scalar, const scalar);

            virtual void boundaryRelaxMatrix
            (
                fvBlockMatrix<vector>&
            ) const override;

            virtual void updateCoeffs() override;


            virtual tmp<Field<vector>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const override;

            virtual tmp<Field<vector>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const override;

            virtual tmp<Field<vector>> valueDivInternalCoeffs
            (
                const tmp<scalarField>&
            ) const override;

            virtual tmp<Field<vector>> valueDivBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const override;

            virtual tmp<Field<vector>> gradientInternalCoeffs() const override;

            virtual tmp<Field<vector>> gradientBoundaryCoeffs() const override;

        //- Write

            virtual void write(Ostream&) const override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
