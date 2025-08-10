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
    (c) 2017-2023 Esi Ltd.

Class
    Foam::flowRateInletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a velocity boundary condition, derived
    from the flux (volumetric or mass-based), whose direction is assumed
    to be normal to the patch.

    For a mass-based flux:
    - the flow rate should be provided in kg/s
    - if \c rho is "none" the flow rate is in m3/s
    - otherwise \c rho should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must
      specify the inlet density using the \c rhoInlet entry

    For a volumetric-based flux:
    - the flow rate is in m3/s

Usage
    \table
        Property     | Description             | Required    | Default value
        massFlowRate | mass flow rate [kg/s]   | no          |
        volumetricFlowRate | volumetric flow rate [m3/s]| no |
        rhoInlet     | inlet density           | no          |
        extrapolateProfile | Extrapolate velocity profile | no | false
        lpmUnits | volumetricFlowRate in liters/min | no | false
    \endtable

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    <patchName>
    {
        type                flowRateInletVelocity;
        volumetricFlowRate  0.2;
        rotating    false;
        axis        (1 0 0);
        origin      (0 0 0);
        angle       constant 10; //angle in degrees between radial vector and velocity
        extrapolateProfile  yes;
        value               uniform (0 0 0);
    }
    \endverbatim

    Example of the boundary condition specification for a mass flow rate:
    \verbatim
    <patchName>
    {
        type                flowRateInletVelocity;
        massFlowRate        0.2;
        extrapolateProfile  yes;
        rho                 rho;
        rhoInlet            1.0;
        rotating    false;
        axis        (1 0 0);
        origin      (0 0 0);
        angle       constant 10; //angle in degrees between radial vector and velocity
        value               uniform (0 0 0);
    }
    \endverbatim

    The \c flowRate entry is a \c Function1 of time, see Foam::Function1Types.

Note
    - \c rhoInlet is required for the case of a mass flow rate, where the
      density field is not available at start-up
    - The value is positive into the domain (as an inlet)
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types

SourceFiles
    flowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateInletVelocityFvPatchVectorField_H
#define flowRateInletVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvPatchFields/derived/flowRateBase/flowRateBase.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class flowRateInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class flowRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public flowRateBase
{
    // Private data

        // Rotational specification

            //- switch rotational component on/off
            Switch rotating_;
            //- axis of rotation
            autoPtr<vector> axis_;
            //- axis origin
            autoPtr<vector> origin_;

        //- resulting flow angle in degrees measured from radial vector between
        //- axis and face centre,
        //- rotational direction uses right hand rule around axis
            autoPtr<Function1<scalar>> angle_;

        //- Set true to extrapolate the velocity profile from the interior
        Switch extrapolateProfile_;

        //- patch relax in block solver
        scalar relax_;

        //- flux correction when using the BC at outlet to make flux value
        // consistent due to the reverse convection coefficients
        Switch fluxCorrection_;


    // Private member functions

public:
        //- Update the patch values given the appropriate density type and value
        template<class RhoType>
        void updateValues(const RhoType& rho);


public:

   //- Runtime type information
   TypeName("flowRateInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        flowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        flowRateInletVelocityFvPatchVectorField
        (
            const flowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        flowRateInletVelocityFvPatchVectorField
        (
            const flowRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const override
        {
            return tmp<fvPatchVectorField>
            (
                new flowRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        flowRateInletVelocityFvPatchVectorField
        (
            const flowRateInletVelocityFvPatchVectorField&,
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
                new flowRateInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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


        //- Return the matrix diagonal coefficients
        virtual tmp<Field<vector>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const override;

        //- Return the matrix source coefficients
        virtual tmp<Field<vector>> valueBoundaryCoeffs
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

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs() override;

        //- Manipulate matrix with given weights
        virtual void boundaryRelaxMatrix(fvBlockMatrix<vector>&) const override;

        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
