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
    (c) 2006 Mark Olesen
    (c) 2017 Esi Ltd.

Class
    fanVelocityFvPatchVectorField

Description
    A fan boundary that derives the mass flow from the pressure change.

    Boundary works as an inlet, outlet or outlet for an explicitly coupled patch
    pair. User should specify either the neighbour pressure or the neighbour
    patch. Fan curve can be specified as volumetric or mass flow based. Switch
    enables user to use Pascal based tables, so no manual conversion is needed.
    filterSize is the number of iterations to average the pressure difference
    over. A larger filter is more stable, but will lagg the solution more.
    The filter is specified via a Function1, so it can vary during the solution.
    Outlet fanVelocity boundaries appear to be more stable than inlet
    boundaries.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        fanCurve     | Volumetric fan curve function | no   |
        //volumetricFanCurve     | Volumetric fan curve function | no   |
        //massFanCurve     | Mass flow fan curve function | no   |

        direction    | direction of flow through fan [in/out] | no | out
        p            | Pressure field name | no | p
        phi          | Flux field name | no | phi
        rho          | Density name | no | rho

        pressureInputType | Neighbour fan pressure [fixed/patch] | no | fixed
        p0           | Fixed neighbour pressure | no    | 0.0
        fanOutletPatch | Patch to obtain neighbour pressure [patch][in] | yes |
        fanInletPatch | Patch to obtain neighbour pressure [patch][out] | yes |

        filterSize    | How many iterations to average the fan pressure drop over | no | 200
        convertPressureToPascal | use pressure in Pa to acces fan curve | no | true

        rhoRef        | Reference density for incompressible cases | yes |
        rhoAverage    | Use patch wise averaging for calculating density | no | false
        maxDelta      | Maximum fractional rate of change of flow rate | no | 0.2

        uniform       | Enforce uniform velocity specification | no | false outlet/true inlet
        fluxValueSwitch | Determines boundary formulation switching based on flux or flow rate | no | true

    \endtable

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            fanVelocity;
        fanCurve        table
        (
            (0     10)
            (100   0)
        );
        outOfBounds     clamp;
        rhoRef          1.2;
        direction       in;
        value           uniform (0 0 0);
    }

    outlet
    {
        type            fanVelocity;
        fanCurve        table
        (
            (0     10)
            (100   0)
        );
        outOfBounds     clamp;
        rhoRef          1.2;
        direction       out;
        filterSize      constant 500;
        p0              uniform 100;
        value           uniform (0 0 0);
    }

    outletPair //requires mappedField velocity on inlet patch
    {
        type            fanVelocity;
        fanCurve        table
        (
            (0     10)
            (100   0)
        );
        outOfBounds     clamp;
        direction       out;
        filterSize      constant 500;
        pressureInputType patch;
        fanInletPatch   fanIn;
        rhoRef          1.2;

        value           uniform (0 0 0);
    }

    \endverbatim


SourceFiles
    fanVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fanVelocityFvPatchVectorField_H
#define fanVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/derived/flowRateOutletVelocity/flowRateOutletVelocityFvPatchVectorField.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "containers/LinkedLists/user/FIFOStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class massFluxInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fanVelocityFvPatchVectorField
:
    public flowRateOutletVelocityFvPatchVectorField
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

        //- Fan pressure input type
        enum fanPressureInputType
        {
            fpiFixed,
            fpiPatch
        };

        //- Fan pressure input type names
        static const NamedEnum<fanPressureInputType, 2>
            fanPressureInputTypeNames_;


    // Private data
private:

        //- Tabulated fan curve
        autoPtr<Function1<scalar>> fanCurve_;

        //- Direction of flow through the fan relative to patch
        fanFlowDirection direction_;

        // pressure field name
        word pName_;

        //- how to determine fan pressure change
        fanPressureInputType pressureInput_;

        // pressure on the other side of the fan
        scalar fanNbPatchPressure_;

        // patch for dynamic pressure calculation
        word pressurePatch_;

        // update indicator
        Switch convertPressureToPascal_;

        //number of pressure deltas to store
        autoPtr<Function1<scalar>> sampleSize_;

        // FIFO stack of historic pressure changes across fan
        FIFOStack<scalar> storedDp_;

        // FIFO stack of historic time step sizes
        FIFOStack<scalar> storedDt_;

       // Make sure it is called only once per iteration
       label curTimeIndex_;

    // Private functions

        // return neighbour patch pressure
        scalar patchMeanPressure(label) const;

        // provide filtered pressure drop from recorded data
        scalar filterInput(scalar val, scalar dt);


public:

    //- Runtime type information
    TypeName("fanVelocity");

    // Constructors

        //- Construct from patch and internal field
        fanVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fanVelocityFvPatchVectorField
        //  onto a new patch
        fanVelocityFvPatchVectorField
        (
            const fanVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fanVelocityFvPatchVectorField
        (
            const fanVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fanVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fanVelocityFvPatchVectorField
        (
            const fanVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (const DimensionedField<vector, volMesh>& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new fanVelocityFvPatchVectorField(*this, iF)
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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
