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
    (c) 2020-2023 Esi Ltd.

Class
    Foam::targetFlowRateOutletPressureFvPatchScalarField

Author
    Georgios Karpouzas

Description
    Foam::targetFlowRateOutletPressureFvPatchScalarField

    Modify pressure to satisfy user-specified target flowrate.
    The relation between pressure-drop and flowrate is approximated based on
    Bernulli's equation. (or kinetic Energy Loss).
    So the pressure is corrected every iteration as follows:

        dp = relax*0.5*rho*((massFlowCur^2 - massFlowTar^2)/(rho*surface^2)

    here

    p: pressure drop
    rho: density
    relax: relaxation factor (0.05 default value)
    massFlowCur: current massFlow
    massFlowTar: target massFlow

Usage
    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    myPatch
    {
        type                      targetFlowRateOutletPressure;
        volumetricFlowRateTarget  constant 1.5;
    }
    \endverbatim

    Example of the boundary condition specification for a mass flow rate:
    \verbatim
    myPatch
    {
        type                targetFlowRateOutletPressure;
        massFlowRateTarget  constant 1.5;
        rho                 rho;
    }
    \endverbatim


SourceFiles
    targetFlowRateOutletPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef targetFlowRateOutletPressureFvPatchScalarField_H
#define targetFlowRateOutletPressureFvPatchScalarField_H

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"
#include "fields/fvPatchFields/derived/flowRateBase/flowRateBase.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class targetFlowRateOutletPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class targetFlowRateOutletPressureFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public flowRateBase
{
    // Private data

        //- Name of the flux
        word phiName_;

        //- Relaxation factor for updating the pressure
        scalar relax_;

        //- Min limiter
        scalar minP_;

        //- Max limiter
        scalar maxP_;

        //- Rho
        scalar rho_;


    // Private member functions

        //- Function checks flux dimensions and returns if
        //  they are massflux
        bool checkPhiDimensions
        (
            const fvsPatchField<scalar>& phiP
        ) const;

        //- Calculate target flow rate
        scalar calcTargetFlowRate
        (
            const bool isPhiMassFlow
        ) const;


public:

    //- Runtime type information
    TypeName("targetFlowRateOutletPressure");


    // Constructors

        //- Construct from patch and internal field
        targetFlowRateOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        targetFlowRateOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  targetFlowRateOutletPressureFvPatchScalarField onto a new patch
        targetFlowRateOutletPressureFvPatchScalarField
        (
            const targetFlowRateOutletPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        targetFlowRateOutletPressureFvPatchScalarField
        (
            const targetFlowRateOutletPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new targetFlowRateOutletPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        targetFlowRateOutletPressureFvPatchScalarField
        (
            const targetFlowRateOutletPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new targetFlowRateOutletPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            //  reversed to true 2022
            virtual bool fixesValue() const
            {
                return true;
            }

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            void updateInletFlowRate(const scalar&);


        //- Write
        virtual void write(Ostream&) const;

        //- Read constant density if incompressible
        scalar calcRho
        (
            const word& rhoName,
            const fvPatchField<scalar>* rhop,
            const scalar surface
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
