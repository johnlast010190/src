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
    (c) 2016 OpenFOAM Foundation
    (c) 2017 Vascular Flow
    (c) 2023 Esi Ltd

Class
    Foam::flowRateTubeVelocityProfileFvPatchVectorField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a velocity boundary condition for circular
    patches, in which a smooth velocity profile is given. The velocity profile
    is scaled by derived the flux (volumetric or mass-based), whose direction is
    assumed to be normal to the patch.

    For a mass-based flux:
    - the flow rate should be provided in kg/s
    - if \c rhoName is "none" the flow rate is in m3/s
    - otherwise \c rhoName should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must
      specify the inlet density using the \c rhoInlet entry

    For a volumetric-based flux:
    - the flow rate is in m3/s

    The velocity profile is given by

        \f[
            U(r) = \frac{Q}{A} \frac{n + 2}{n}
            \left[ 1 - \left(\frac{r}{R} \right)^n \right]
        \f]

    where
    \vartable
        U      | velocity normal to patch [m/s]
        Q      | volumetricflow rate [m^3/s]
        A      | circular patch area [m^2]
        R      | circular patch radius [m]
        r      | radial direction from patch centre [m]
        n      | degree of polynomial profile [-]
    \endtable

    \heading Patch usage

    \table
        Property   | Description                    | Req'd?    | Default
        massFlowRate | mass flow rate [kg/s]        | no        |
        volumetricFlowRate  | volumetric flow rate  | no        |
        rhoInlet     | inlet density                | no        |
        n          | degree of polynomial profile   | no        | 2
    \endtable

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    myPatch
    {
        type    flowRateTubeVelocityProfile;
        n       6;
        volumetricFlowRate table
        (
            (0 1)
            (1 2)
            (2 1)
        );
        outOfBounds clamp;
    }
    \endverbatim

    Example of the boundary condition specification for a mass flow rate:
    \verbatim
    myPatch
    {
        type                flowRateInletVelocity;
        massFlowRate        0.2;
        rho                 rho;
        rhoInlet            1.0;
    }
    \endverbatim

    The \c flowRate entry is a \c Function1 type, meaning that it can be
    specified as constant, a polynomial fuction of time, and ...

Note
    - The patch geometry must be circular.
    - Parabolic profile is given by n = 2 (default).
    - \c rhoInlet is required for the case of a mass flow rate, where the
      density field is not available at start-up
    - the value is positive into the domain (as an inlet)
    - may not work correctly for transonic inlets
    - strange behaviour with potentialFoam since the U equation is not solved

SeeAlso
    Foam::Function1Types
    Foam::fixedValueFvPatchField

SourceFiles
    flowRateTubeVelocityProfileFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateTubeVelocityProfileFvPatchVectorField_H
#define flowRateTubeVelocityProfileFvPatchVectorField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvPatchFields/derived/flowRateBase/flowRateBase.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class flowRateTubeVelocityProfileFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class flowRateTubeVelocityProfileFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public flowRateBase
{
    // Private data

        //- degree of polynomial profile
        scalar n_;


public:

    //- Runtime type information
    TypeName("flowRateTubeVelocityProfile");


    // Constructors

        //- Construct from patch and internal field
        flowRateTubeVelocityProfileFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flowRateTubeVelocityProfileFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        flowRateTubeVelocityProfileFvPatchVectorField
        (
            const flowRateTubeVelocityProfileFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        flowRateTubeVelocityProfileFvPatchVectorField
        (
            const flowRateTubeVelocityProfileFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new flowRateTubeVelocityProfileFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        flowRateTubeVelocityProfileFvPatchVectorField
        (
            const flowRateTubeVelocityProfileFvPatchVectorField&,
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
                new flowRateTubeVelocityProfileFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<vector>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<vector>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

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
