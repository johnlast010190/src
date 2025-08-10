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
    (c) 2010-2023 Esi Ltd.

Class
    Foam::flowRateOutletVelocityFvPatchVectorField

Description
    An outlet boundary for velocity that scales the normal component of
    a zero gradient boundary such that the overall massflow matches the
    specified value

    Includes a "uniform" option which switches to a fixed uniform formulation
    with a "weak" convection boundary implementation

    For a mass-based flux:
    - the flow rate should be provided in kg/s
    - if \c rhoName is "none" the flow rate is in m3/s
    - otherwise \c rhoName should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must
      specify the inlet density using the \c rhoInlet entry

    For a volumetric-based flux:
    - the flow rate is in m3/s


    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    myPatch
    {
        type        flowRateOutletVelocity;
        volumetricFlowRate  0.2;
        value       uniform (0 0 0); // placeholder
    }
    \endverbatim

    Example of the boundary condition specification for a mass flow rate:
    \verbatim
    myPatch
    {
        type                flowRateOutletVelocity;
        massFlowRate        0.2;
        rho                 rho;
        rhoInlet            1.0;
    }
    \endverbatim


SourceFiles
    flowRateOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateOutletVelocityFvPatchVectorFields_H
#define flowRateOutletVelocityFvPatchVectorFields_H

#include "fields/fvPatchFields/derived/fixedNormalSlip/fixedNormalSlipFvPatchFields.H"
#include "fields/fvPatchFields/derived/flowRateBase/flowRateBase.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class flowRateOutletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class flowRateOutletVelocityFvPatchVectorField
:
    public fixedNormalSlipFvPatchVectorField,
    public flowRateBase
{
    // Private data

        //- Outlet integral flow rate
        autoPtr<Function1<scalar>> flowRate_;

protected:
        //- current flowrate for switches, positive=outlet, negative=inlet
        scalar currentFlowRate_;

        //- maximum rate of change
        scalar maxDelta_;

        //- flux name
        word phiName_;

        // switch - use uniform outlet velocity with weak convection BC
        Switch uniform_;

        // switch for determining switching mechanism for boundary direction
        Switch fluxValueSwitch_;

        //- patch relax in block solver
        scalar relax_;

private:
    // Private functions

        //- scale provided flux to match flowRate
        void scaleFlux(const scalar, scalarField&) const;

        //- calculate non-uniform velocity at outlet patch
        tmp<scalarField> nonUniformVelocity(const scalar) const;

        //- calculate uniform velocity at inlet or outlet
        tmp<scalarField> uniformVelocity(const scalar) const;

protected:
    // Protected functions

        //- dictionary constructor checks
        void initialise(const dictionary& dict);

public:

    //- Runtime type information
    TypeName("flowRateOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        flowRateOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flowRateOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateOutletVelocityFvPatchVectorField onto a new patch
        flowRateOutletVelocityFvPatchVectorField
        (
            const flowRateOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        flowRateOutletVelocityFvPatchVectorField
        (
            const flowRateOutletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const override
        {
            return tmp<fvPatchVectorField>
            (
                new flowRateOutletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        flowRateOutletVelocityFvPatchVectorField
        (
            const flowRateOutletVelocityFvPatchVectorField&,
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
                new flowRateOutletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<vector>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const override;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
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

        //- Update the coefficients associated with the patch field given
        //  uniform density field
        void updateCoeffs(const scalar uniformRho);

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs() override;

        //- Manipulate matrix with given weights
        virtual void boundaryRelaxMatrix
        (
            fvBlockMatrix<vector>&
        ) const override;

        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
