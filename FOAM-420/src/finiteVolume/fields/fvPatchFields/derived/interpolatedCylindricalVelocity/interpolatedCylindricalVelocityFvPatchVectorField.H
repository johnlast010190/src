/*---------------------------------------------------------------------------*\
| Modified 2010-2012 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
                    |
     F ield         | FOAM: The Open Source CFD Toolbox
     O peration     |
     A nd           | Copyright (C) 2009 Icon CG Ltd.
     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of FOAM.

    FOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    FOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::interpolatedCylindricalVelocityFvPatchVectorField

Description
    A table based radialy-varying velocity boundary condition.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            interpolatedCylindricalVelocity;
        phi             phi;
        rho             rho;

        profile table   //use Function1 interface
        (
            (r1 (Uz1 Ur1 Ut1))
            (r2 (Uz2 Ur2 Ut2))
            ...

            (rN (UzN UrN UtN))
        );
        outOfBounds     clamp;    // (error|warn|clamp|repeat)

        centre           0 0 0); //centre of cylindrical coordinates
        axis            (0 0 1);

        flowRate        1; //kg/s for compressible and m3/s for incompressible
                         // this entry is optional and used to scale normal
                         // velocity distribution if present
    }
    @endverbatim

See Also
    Foam::interpolationTable and Foam::fixedValueFvPatchField

SourceFiles
    interpolatedCylindricalVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef interpolatedCylindricalVelocityFvPatchVectorField_H
#define interpolatedCylindricalVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchField.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class interpolatedCylindricalVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class interpolatedCylindricalVelocityFvPatchVectorField
:
    public fixedValueFvPatchField<vector>
{

private:
    // Private data

        //- The distribution being used, including the bounding treatment
        autoPtr<Function1<vector>> distribution_;

        //- flux name, "none" for fixedValue
        word phiName_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        //- point on cylindrical coordinate axis
        point centre_;

        //- cylinder axis
        vector axis_;

        //- optional flow rate specification
        autoPtr<scalar> flowRate_;

        //- Private member functions

            //- assign values from interpolation table
            void mapTable();

public:

    //- Runtime type information
    TypeName("interpolatedCylindricalVelocity");


    // Constructors

        //- Construct from patch and internal field
        interpolatedCylindricalVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interpolatedCylindricalVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        interpolatedCylindricalVelocityFvPatchVectorField
        (
            const interpolatedCylindricalVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        interpolatedCylindricalVelocityFvPatchVectorField
        (
            const interpolatedCylindricalVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return tmp<fvPatchField<vector>>
            (
                new interpolatedCylindricalVelocityFvPatchVectorField
                (*this)
            );
        }

        //- Construct as copy setting internal field reference
        interpolatedCylindricalVelocityFvPatchVectorField
        (
            const interpolatedCylindricalVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector>>
            (
                new interpolatedCylindricalVelocityFvPatchVectorField
                (*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
