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
    (c) 2017 OpenCFD Ltd
    (c) 2022 Esi Ltd.

Class
    Foam::fanPressureFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition can be applied to assign either a static/total
    pressure inlet or outlet for a fan. In case of a patch blowing out
    with total pressure specification the boundary condition is
    switched to static pressure based on flux (dynamic component is not added).

    The switch nonDimensional can be used for a non-dimensional table. It needs
    inputs rpm and dm of the fan.

    The nonDimensional flux for the table is calculate as :

        phi = 4.0*mDot/(rho*sqr(PI)*dm^3*omega)
        where:
            dm is the mean diameter.
            omega is rad/sec.

    The nonDimensinal pressure :

        Psi = 2 deltaP/(rho*(sqr(PI*omega*dm)))
        where:
            deltaP is the pressure drop

    The non-dimensional table should be given as Psi = F(phi).

Usage
    \table
        Property     | Description             | Required    | Default value
        file         | fan curve file name     | no         |
        outOfBounds  | out of bounds handling  | no         |
        direction    | direction of flow through fan [in/out] | yes |
        pressureMode | static or total pressure [static/total] | no | static
        p0           | environmental total pressure | yes    |
        nonDimensional | uses non-dimensional table | no | false
        rpm          | fan rpm for non-dimensional table     | no | 0.0
        dm           | mean diameter for non-dimensional table | no | 0.0
        timeFanCurve | is table time depenent? | no          | false
    \endtable

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            fanPressure;
        file            "fanCurve";
        outOfBounds     clamp;
        direction       in;
        pressureMode    total;
        p0              uniform 0;
        value           uniform 0;
        timeFanCurve    false;
    }

    outlet
    {
        type            fanPressure;
        file            "fanCurve";
        outOfBounds     clamp;
        direction       out;
        p0              uniform 0;
        value           uniform 0;
        timeFanCurve    false;
    }

    inlet
    {
        type            fanPressure;
        timeFanCurve    true;
        type            uniformTable;
        low             (0  0);
        high            (10 30);
        values 3 4
        (
            (11 11 12 12)
            (13 13 14 14)
            (15 15 16 16)
        );
        direction       in;
        pressureMode    static;
        p0              uniform 0;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::fanFvPatchField
    Foam::totalPressureFvPatchScalarField
    Foam::interpolationTable

SourceFiles
   fanPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fanPressureFvPatchScalarField_H
#define fanPressureFvPatchScalarField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/functions/Function2/Function2/Function2.H"
#include "addStaticHead/staticHead.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fanPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fanPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

public:

        //- Fan flow direction
        enum fanFlowDirection
        {
            ffdIn,
            ffdOut
        };

        //- Fan flow directions names
        static const NamedEnum<fanFlowDirection, 2> fanFlowDirectionNames_;

        //- pressure modes static/total
        enum pressureMode
        {
            pSt,
            pTot
        };

        //- pressure mode names
        static const NamedEnum<pressureMode, 2> pressureModeNames_;


private:

    // Private data

        //- Old value of the patch field
        autoPtr<scalarField> oldPressure_;

        //- Stored time index
        label timeIndex_;

        //- Update ones per timestep
        Switch updateOnesPerTimestep_;

        //- Relaxation factor
        scalar relax_;

        //- Is the fan curve dependent on time
        Switch isTimeCurve_;

        //- Tabulated fan curve
        autoPtr<Function1<scalar>> fanCurve_;

        //- Tabulated dependent on time/flow rate
        autoPtr<Function2<scalar>> timeFanCurve_;

        //- Direction of flow through the fan relative to patch
        fanFlowDirection direction_;

        //- Swtich for using non-dimensional curve
        Switch nonDimensional_;

        // Parameters for non-dimensional table

            //- Fan rpm
            scalar rpm_;

            //- Fan mean diameter
            scalar dm_;

         //- Total pressure parameters
            pressureMode pressureMode_;

            //- Name of the velocity field
            word UName_;

            //- Name of the flux transporting the field
            word phiName_;

            //- Total pressure
            scalarField p0_;

            //- Static head calculator type class
            staticHead staticHead_;

public:

    //- Runtime type information
    TypeName("fanPressure");


    // Constructors

        //- Construct from patch and internal field
        fanPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fanPressureFvPatchScalarField
        //  onto a new patch
        fanPressureFvPatchScalarField
        (
            const fanPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fanPressureFvPatchScalarField
        (
            const fanPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fanPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fanPressureFvPatchScalarField
        (
            const fanPressureFvPatchScalarField&,
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
                new fanPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Inherit updateCoeffs from fixedValueFvPatchScalarField
        using fixedValueFvPatchScalarField::updateCoeffs;

        //- Add dynamic (if needed) and static head
        virtual void addDynamicAndStaticHead
        (
            const scalarField& p0p,
            const vectorField& Up
        );

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
