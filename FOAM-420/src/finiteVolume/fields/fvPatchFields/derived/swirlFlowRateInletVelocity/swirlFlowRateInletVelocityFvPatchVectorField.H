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
    (c) 2023 Esi Ltd.

Class
    Foam::swirlFlowRateInletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a volumetric- OR mass-flow normal vector
    boundary condition by its magnitude as an integral over its area with a
    swirl component determined by the angular speed, given in revolutions per
    minute (RPM)

    The basis of the patch (volumetric or mass) is determined by the
    dimensions of the flux, phi. The current density is used to correct the
    velocity when applying the mass basis.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        flowRate     | flow rate profile       | yes         |
        rpm          | rotational speed profile | yes        |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            swirlFlowRateInletVelocity;
        flowRate        constant 0.2;
        rpm             constant 100;
    }
    \endverbatim

Note
    - the \c flowRate and \c rpm entries are Function1 types, able to describe
      time varying functions.  The example above gives the usage for supplying
      constant values.
    - the value is positive into the domain

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types

SourceFiles
    swirlFlowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef swirlFlowRateInletVelocityFvPatchVectorField_H
#define swirlFlowRateInletVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvPatchFields/derived/flowRateBase/flowRateBase.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
        Class swirlFlowRateInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class swirlFlowRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public flowRateBase
{
    // Private data

        //- Angular speed in revolutions per minute (RPM)
        autoPtr<Function1<scalar>> rpm_;

        // enable weak formulation for outlets
        Switch outlet_;

public:

   //- Runtime type information
   TypeName("swirlFlowRateInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const swirlFlowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const swirlFlowRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new swirlFlowRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const swirlFlowRateInletVelocityFvPatchVectorField&,
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
                new swirlFlowRateInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<vectorField> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<vectorField> valueBoundaryCoeffs
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
