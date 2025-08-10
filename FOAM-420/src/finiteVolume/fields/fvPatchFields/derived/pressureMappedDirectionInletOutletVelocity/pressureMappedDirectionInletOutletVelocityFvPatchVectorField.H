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
    (c) 2020 Esi Ltd.

Class
    Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    A pressureDirectedInletOutletVelocity with mapped direction of the velocity
    components.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        mapMethod    | Type of mapping         | no       | planarInterpolation
        perturb      | Perturb points for regular geometries | no | 1e-5
        delabella      | Whether to use delabella triangulation | no | true
        cylindricalCoords | Switch to convert values from cylindrical to
            cartesian coordinate system | no | false
        preserveSign    | preserve the sign of directrion | no       | false
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            pressureMappedDirectionInletOutletVelocity;
        phi             phi;
        rho             rho;
        perturb         0;
        delabella         true;
        mapMethod       nearest;
        value           uniform (0 0 0);
    }
    \endverbatim

Note
    Sign conventions:
    - positive flux (out of domain): apply zero-gradient condition
    - negative flux (into of domain): derive from the flux with specified
      direction

See also
    Foam::mixedFvPatchVectorField

SourceFiles
    pressureMappedDirectionInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureMappedDirectionInletOutletVelocityFvPatchVectorField_H
#define pressureMappedDirectionInletOutletVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"
#include "triSurface/triSurfaceTools/pointToPointPlanarInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class pressureMappedDirectionInletOutletVelocity Declaration
\*---------------------------------------------------------------------------*/

class pressureMappedDirectionInletOutletVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private data

        //- Name of the field data table, defaults to the name of the field
        word fieldTableName_;

        //- Flux field name
        word phiName_;

        //- Density field name
        word rhoName_;

        //- Inlet direction
        vectorField inletDir_;

        //- Fraction of perturbation (fraction of bounding box) to add
        scalar perturb_;

        //- Whether to use delabella triangulation
        Switch delabella_;

        //- Interpolation scheme to use
        word mapMethod_;

        //- If true, points are in cylindrical coordinates and must be
        //  converted to cartesian
        Switch cylindricalCoords_;

        //- 2D interpolation (for 'planarInterpolation' mapMethod)
        autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- If true, force velocity to be at the same side with inletDir_
        Switch preserveSign_;

    // Private Member functions

        //- Check and do mapping
        void checkAndMap();

public:

    //- Runtime type information
    TypeName("pressureMappedDirectionInletOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        pressureMappedDirectionInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureMappedDirectionInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  pressureMappedDirectionInletOutletVelocityFvPatchVectorField
        //  onto a new patch
        pressureMappedDirectionInletOutletVelocityFvPatchVectorField
        (
            const pressureMappedDirectionInletOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureMappedDirectionInletOutletVelocityFvPatchVectorField
        (
            const pressureMappedDirectionInletOutletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new pressureMappedDirectionInletOutletVelocityFvPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        pressureMappedDirectionInletOutletVelocityFvPatchVectorField
        (
            const pressureMappedDirectionInletOutletVelocityFvPatchVectorField&,
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
                new pressureMappedDirectionInletOutletVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Attributes

            //- Return true: this patch field is altered by assignment
            virtual bool assignable() const
            {
                return true;
            }


        // Access

            //- Return the name of rho
            const word& rhoName() const
            {
                return rhoName_;
            }

            //- Return reference to the name of rho to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
            }


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


        // Utility functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
