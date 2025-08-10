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
    (c) 2019 Esi Ltd.

Class
    Foam::wallVelocityFvPatchVectorField

Group

Description
    Wall boundary that inforces zero wall normal normal stress.

SourceFiles
    wallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef wallVelocityFvPatchVectorField_H
#define wallVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class wallVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class wallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("wallVelocity");


    // Constructors

        //- Construct from patch and internal field
        wallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wallVelocityFvPatchVectorField
        //  onto a new patch
        wallVelocityFvPatchVectorField
        (
            const wallVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wallVelocityFvPatchVectorField
        (
            const wallVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const override
        {
            return tmp<fvPatchVectorField>
            (
                new wallVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wallVelocityFvPatchVectorField
        (
            const wallVelocityFvPatchVectorField&,
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
                new wallVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const override
            {
                return true;
            }

            //- Is defined in frame supported?
            virtual bool definedInFrameSupported() const override
            {
                return true;
            }

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs() override;

        // Evaluation functions

            virtual tmp<vectorField> gradientBoundaryCoeffs() const override;

            virtual tmp<CoeffField<vector>>
            gradientInternalBCoeffs() const override;

            virtual tmp<vectorField> gradientBoundaryBCoeffs() const override;

            //- Return gradient at boundary
            virtual tmp<vectorField> snGrad() const override;

        //- Write
        virtual void write(Ostream&) const override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
