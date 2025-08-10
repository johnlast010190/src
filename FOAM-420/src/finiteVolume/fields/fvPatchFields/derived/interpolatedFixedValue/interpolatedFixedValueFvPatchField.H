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
    Foam::interpolatedFixedValueFvPatchField

Description
    A spatially-varying form of a fixed value boundary condition.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            interpolatedFixedValue;
        fileName        "$FOAM_CASE/spatial-distribution"; // optional required data
        data
        {
            x1 y1
            x2 y2
            ...

            xN yN
        }               //local data used if fileName entry absent
        outOfBounds     clamp;           // (error|warn|clamp|repeat)
        distanceType    z; // (x, y, z, wallDistance (distance from wall),
                        // (cont.) distanceToPoint, distanceAlongVector
        point           (0  0 0); //entry required for pointDistance and
                        // distanceAlongVector
        direction       (1 0 0); //entry required for distanceAlongVector
        runTimeModifiable false; // default false - determines whether
                                    the data file
                        // is re-read after each iteration
                        // allows limited (and unreliable) external coupling
        xcale           1;
        yscale          1;
        xoffset         0;
        yoffset         Type::zero;
        //y = yoffset + yscale * F((x - xoffset)/xscale)
    }
    @endverbatim

Note
    This class is derived directly from a fixedValue patch rather than from
    a uniformFixedValue patch.

See Also
    Foam::Function1 and Foam::fixedValueFvPatchField

SourceFiles
    interpolatedFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef interpolatedFixedValueFvPatchField_H
#define interpolatedFixedValueFvPatchField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchField.H"
#include "cfdTools/general/fieldProfile/fieldProfile.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class interpolatedFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolatedFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{

private:
    // Private data

        //- The distribution being used, including the bounding treatment
        autoPtr<fieldProfile<Type>> distribution_;

public:

    //- Runtime type information
    TypeName("interpolatedFixedValue");


    // Constructors

        //- Construct from patch and internal field
        interpolatedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interpolatedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        interpolatedFixedValueFvPatchField
        (
            const interpolatedFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        interpolatedFixedValueFvPatchField
        (
            const interpolatedFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const override
        {
            return tmp<fvPatchField<Type>>
            (
                new interpolatedFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        interpolatedFixedValueFvPatchField
        (
            const interpolatedFixedValueFvPatchField<Type>&,
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
                new interpolatedFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs() override;

            //- Approximate dY/dZ for face centres
            tmp<Field<Type>> gradient() const;

            //- Is defined in frame supported?
            virtual bool definedInFrameSupported() const override
            {
                return true;
            }

        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
