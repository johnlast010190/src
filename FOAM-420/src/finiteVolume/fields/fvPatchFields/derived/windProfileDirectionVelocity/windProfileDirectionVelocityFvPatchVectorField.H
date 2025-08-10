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
    (c) 2010-2014 Esi Ltd.
    (c) 1991-2009 OpenCFD Ltd.

Class
    Foam::windProfileDirectionVelocityFvPatchVectorField

Description
Description
    This boundary condition provides spatially 1D varying inletOutlet
    that can be rotated around a user specified axis in a time dependent manner.
    Inflow velocity magnitude profile is specified via a Function1: profile.
    Input-coordinate is supplied via patchDistanceFunction.

    Direction of wind is specified in degrees around an axis.
    The direction vector is calculated around Z and then rotated into
    the specified axis frame. The base Z axis frame has the following
    physical relationship to wind direction:

        0deg = -Y = (0 -1 0) = North
        90deg = -X = (-1 0 0) =  East
        180deg = Y = (0 1 0) =  South
        270deg = X = (1 0 0) = West


    The boundary also provides velocity gradient data for turbulent dissipation
    calculation

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        profile      | wind velocity magnitude | yes         |
        xscale       | scaling for distance    | no          | 1
        yscale       | scaling for profile     | no          | 1
        xoffset      | distance additive const | no          | 0
        yoffset      | profile additive        | no          | 0
        distanceType | distance function       | yes         |
        windDirection| wind direction [degrees]| yes         |
        axis         | direction axis          | no          | (0 0 1)
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            windProfileDirectionVelocity;
        profile         table
        (
            (0  1)
            (1  1.5)
            (5  2)
        );
        distanceType    wallDistance;
        axis            (0 0 1);
        windDirection   constant 45;
        normalDirection false;

        value uniform   (1.5 0 0);
    }
    \endverbatim


SourceFiles
    windProfileDirectionVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef windProfileDirectionVelocityFvPatchVectorField_H
#define windProfileDirectionVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"
#include "cfdTools/general/fieldProfile/fieldProfile.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class windProfileDirectionVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class windProfileDirectionVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
private:
    // Private data

        //- The distribution being used, including the bounding treatment
        autoPtr<fieldProfile<scalar>> distribution_;

        //- Flux name, "none" for fixedValue
        word phiName_;

        //- Patch normal direction
        Switch normalDirection_;

        //- Wind direction
        autoPtr<Function1<scalar>> windDirectionDeg_;


    // Private member functions

        vector direction() const;

public:

    //- Runtime type information
    TypeName("windProfileDirectionVelocity");


    // Constructors

        //- Construct from patch and internal field
        windProfileDirectionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        windProfileDirectionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  windProfileDirectionVelocityFvPatchVectorField
        //  onto a new patch
        windProfileDirectionVelocityFvPatchVectorField
        (
            const windProfileDirectionVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        windProfileDirectionVelocityFvPatchVectorField
        (
            const windProfileDirectionVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const override
        {
            return tmp<fvPatchVectorField>
            (
                new windProfileDirectionVelocityFvPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        windProfileDirectionVelocityFvPatchVectorField
        (
            const windProfileDirectionVelocityFvPatchVectorField&,
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
                new windProfileDirectionVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs() override;

        //- Calculate dY/dX for face centres
        tmp<Field<scalar>> gradient() const;

        //- Is defined in frame supported?
        virtual bool definedInFrameSupported() const override
        {
            return true;
        }

        //- Write
        virtual void write(Ostream&) const override;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf) override;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
