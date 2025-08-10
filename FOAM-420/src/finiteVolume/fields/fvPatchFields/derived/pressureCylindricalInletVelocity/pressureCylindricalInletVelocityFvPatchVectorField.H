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
    (c) Esi Ltd
    (c) 2011-2012 OpenFOAM Foundation

Class
    Foam::pressureCylindricalInletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This velocity inlet boundary condition is applied to patches where the
    pressure is specified.  The inflow velocity is obtained from the flux with
    the specified cylindrical coordinate inlet direction.

    \heading Patch usage

    \table
        Property     | Description                    | Required | Default value
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        origin       | cylindrical coordinate origin  | no       | (0 0 0)
        axis         | cylindrical coordinate Z axis  | no       | (0 0 1)
        cylindricalInletDirection | inlet direction in cylindrical coordinate (r, theta[deg], z) per patch face | yes |
        relax        | tangential velocity update relaxation | no | 0.25
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type                        pressureCylindricalInletVelocity;
        phi                         phi;
        rho                         rho;
        cylindricalInletDirection   uniform (36 -60 0);
        relax                       0.1;
        value                       uniform (0 0 0);
    }
    \endverbatim

Note
    Reverse flow is not supported by this boundary type.
    Relaxation mainly required for compressible flow at large angles.


SeeAlso
    Foam::pressureDirectedInletVelocityFvPatchVectorField
    Foam::pressureDirectedInletOutletVelocityFvPatchVectorField

SourceFiles
    pressureCylindricalInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureCylindricalInletVelocityFvPatchVectorField_H
#define pressureCylindricalInletVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/derived/pressureDirectedInletVelocity/pressureDirectedInletVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class pressureCylindricalInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class pressureCylindricalInletVelocityFvPatchVectorField
:
    public pressureDirectedInletVelocityFvPatchVectorField
{
    // Private data

        //- Central point
        const vector origin_;

        //- Axis
        const vector axis_;

        //- Cylindrical vector component
        vectorField cylindricalVector_;

        //- Tangential inlet velocity component update relaxation factor
        scalar relax_;

    // Private member functions

        void convertCylindricalToCartesian();

public:

    //- Runtime type information
    TypeName("pressureCylindricalInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        pressureCylindricalInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureCylindricalInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  pressureCylindricalInletVelocityFvPatchVectorField
        //  onto a new patch
        pressureCylindricalInletVelocityFvPatchVectorField
        (
            const pressureCylindricalInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureCylindricalInletVelocityFvPatchVectorField
        (
            const pressureCylindricalInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new pressureCylindricalInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureCylindricalInletVelocityFvPatchVectorField
        (
            const pressureCylindricalInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new pressureCylindricalInletVelocityFvPatchVectorField
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
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);


        //- Update
        virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
