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
    (c) 2010-2021 Esi Ltd.

Class
    lumpedServerOutletTemperatureFvPatchScalarField

Description
    Compact server boundary condition. Lumped model for outlet
    temperature.
    Outlet patch temperature set as Inlet temperature + offset.

    "A Compact Server Model for Transient Data Center Simulations"
    JamesW. VanGilder, PE Zachary M. Pardey, Member ASHRAE,
    Christopher M. Healey, PhD Xuanhang Zhang
    2013 ASHRAE

    Example of the boundary condition specification:
    @verbatim
    patch
    {
        T
        {
            type lumpedServerOutletTemperature;
            value uniform 300.0;
            patches ( s_in ); //Server inlet patch nbr
            mass 15;  // Server Mass [kg]
            Cps 460; //  Heat Capacity [J/kgK]
            effectiveness 0.8; // Effectiveness [-]
            lambda 0.1; // Geometric factor [-]
            power constant 150; server Power [W]
            massFlowRate 0.015; // Initial value (updated) [kg/s]
            verbose false;
        }
    }
    @endverbatim


SourceFiles
    lumpedServerOutletTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedServerOutletTemperatureFvPatchScalarField_H
#define lumpedServerOutletTemperatureFvPatchScalarField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
         lumpedServerOutletTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class lumpedServerOutletTemperatureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data
private:

        List<wordRe> sourcePatchRegExp_;

        //- Corresponding patchIDs
        labelHashSet patchIDs_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Flux name
        word phiName_;

        //- Additional output
        bool verbose_;

        //- Server data
        //- Mass [kg]
        scalar Ms_;

        //- Specific Heat [J/kgK]
        scalar Cps_;

        //- Effectiveness [-]
        scalar effectiveness_;

        //- Geometric factor [-]
        scalar lambda_;

        //- Power [W]
        autoPtr<Function1<scalar>> Ps_;

        //- Cooling mass flow rate [kg/s]
        scalar qs_;

        //- Server temperature at previous time (old)
        scalar Teff_;

        //- Offset [K]
        scalar deltaT_;

    // Private functions

        //- Return server temperature
        scalar calculateTeff
        (
            const  scalar dt,
            const  scalar tau1,
            const  scalar Tamb,
            const  scalar dTit
        );

        //- Return server delta temperature
        scalar calculateDeltaT
        (
            const  scalar dt,
            const  scalar tau1,
            const  scalar tau2,
            const  scalar Tamb,
            const  scalar dTit
        );

public:

    //- Runtime type information
    TypeName("compressible::lumpedServerOutletTemperature");

    // Constructors

        //- Construct from patch and internal field
        lumpedServerOutletTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        lumpedServerOutletTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  lumpedServerOutletTemperatureFvPatchScalarField
        //  onto a new patch
        lumpedServerOutletTemperatureFvPatchScalarField
        (
            const lumpedServerOutletTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        lumpedServerOutletTemperatureFvPatchScalarField
        (
            const lumpedServerOutletTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new lumpedServerOutletTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        lumpedServerOutletTemperatureFvPatchScalarField
        (
            const lumpedServerOutletTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (const DimensionedField<scalar, volMesh>& iF) const
        {
            return tmp<fvPatchScalarField>
            (
                new lumpedServerOutletTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
