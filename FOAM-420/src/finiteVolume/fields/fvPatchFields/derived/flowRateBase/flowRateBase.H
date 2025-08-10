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
    \endtable

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

#ifndef flowRateBase_H
#define flowRateBase_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "primitives/enums/Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class flowRateBase Declaration
\*---------------------------------------------------------------------------*/

class flowRateBase
{
    // Public data
    public:

        enum massFlowUnits
        {
            kg_s,
            kg_min
        };

        enum volumetricFlowUnits
        {
            m3_s,
            m3_hr,
            l_min,
            l_s
        };

    // Protected data
    protected:

        static const Enum<massFlowUnits> massFlowUnitNames_;
        static const Enum<volumetricFlowUnits> volumetricFlowUnitNames_;

    // Private data
    private:

        //- Accessing the base object registry of the mesh and its fields
        const Foam::objectRegistry& db_;

        //- Reference to the patch class that instantiated flowRateBase
        const fvPatch& patch_;

        //- Is volumetric?
        bool volumetric_;

        //- Inlet integral flow rate
        autoPtr<Function1<scalar>> flowRate_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        //- Rho initialisation value (for start; if value not supplied)
        scalar rhoInlet_;

        //- switch - use average rho to prevent feedback spike
        Switch rhoAverage_;

        //- Flow rate units
        massFlowUnits massFlowUnits_;
        volumetricFlowUnits volumetricFlowUnits_;

        //- Patches selected for inflow when calculating mass flow split
        labelHashSet inflowPatchSet_;

        //- Was the patch list present in the input?
        bool inflowPatchesSpecified_;

        //- Current inlet flux. Used only in split BC mode
        scalar curInletMassFlow_;

        //- for multiphase
        word phaseName_;

        //- Use '...Target' keywords (for pressure boundaries)
        bool targetType_;


    // Constructors

    public:

        TypeName("flowRateBase");

        //- Null constructor
        flowRateBase
        (
            const fvPatch&,
            const Foam::objectRegistry&,
            const bool targetType
        );

        //- Construct from components
        flowRateBase
        (
            const fvPatch&,
            const Foam::objectRegistry&,
            const dictionary&,
            const bool targetType
        );

        //- Copy construct
        flowRateBase
        (
            const fvPatch&,
            const Foam::objectRegistry&,
            const flowRateBase&,
            const bool targetType
        );


    // Destructor

       virtual ~flowRateBase();


    // Public member functions

        //- Calculate current mass flowrate (if !volumetric()) or volumetric
        //  flow rate in SI units, i.e. kg/s or m^3/s
        scalar currentFlowRate() const;

        scalar averageRhoOrOne() const;

        tmp<scalarField> multiplyByRhoOrOne
        (
            const scalarField&, bool useRhoAverage=true
        ) const;

        //- Called to set the inflow rate based on inflow at patches
        void updateInletFlowRate(const List<scalar>& inletFluxes);

        word flowRateKeyword(const word& name) const
        {
            if (targetType_)
            {
                return name + "Target";
            }
            else
            {
                return name;
            }
        }

    // Access functions

        bool volumetric() const
        {
            return volumetric_;
        }

        bool isSplitBC() const
        {
            return (flowRate_->name() == flowRateKeyword("massFlowSplit"));
        }

        word rhoName() const
        {
            return rhoName_;
        }

        Switch rhoAverage() const
        {
            return rhoAverage_;
        }

        scalar rhoInlet() const
        {
            return rhoInlet_;
        }

        word phaseName() const
        {
            return phaseName_;
        }


    //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //