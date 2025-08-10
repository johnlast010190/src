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
    (c) 2010-2012 Esi Ltd.
    (c) 2006 Mark Olesen

Class
    wheelVelocityFvPatchVectorField

Description
    A rotating boundary for a wheel based on the translational velocity of the
    hub and the contact radius.

    Example of the boundary condition specification:
    @verbatim
    wheel
    {
        type            wheelVelocity;
        origin          (0 0 0);
        axis            (0 -1 0);
        hubSpeed        10;
        contactRadius   1;
        value   uniform (0 0 0);
    }
    @endverbatim


SourceFiles
    wheelVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef wheelVelocityFvPatchVectorField_H
#define wheelVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class massFluxInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class wheelVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- axis of rotation
        vector origin_;
        vector axis_;
        scalar hubSpeed_;
        scalar contactRadius_;

    // Private functions

        //- Calculate patch velocity
        void calcVelocity();


public:

    //- Runtime type information
    TypeName("wheelVelocity");

    // Constructors

        //- Construct from patch and internal field
        wheelVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wheelVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  wheelVelocityFvPatchVectorField
        //  onto a new patch
        wheelVelocityFvPatchVectorField
        (
            const wheelVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wheelVelocityFvPatchVectorField
        (
            const wheelVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new wheelVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wheelVelocityFvPatchVectorField
        (
            const wheelVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (const DimensionedField<vector, volMesh>& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new wheelVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

        // Return the origin value
        vector origin() const
        {
            return origin_;
        };

        // Return the axis (of rotation)
        vector axis() const
        {
            return axis_;
        };

        // Return the hubspeed value
        scalar hubSpeed() const
        {
            return hubSpeed_;
        };

        // Return the contactRadius value
        scalar contactRadius() const
        {
            return contactRadius_;
        };

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
