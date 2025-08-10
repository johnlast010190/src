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
    (c) 2022 Esi Ltd.

Class
    Foam::nutkDemirelRoughWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions for rough walls, based on turbulence kinetic
    energy.  The condition manipulates the E parameter to account for roughness
    effects.

    Parameter ranges
    - Roughness height = Equivalent sand grain roughness height (0-10000)

    \table
        Description of condition            | NSTM rating  |  roughness height (micro meter)    |  Rt50 (micro meter)
        Hydraulically smooth surface        |   0          |  0     |   0
        Typical as applied AF coating       |   0          |  30    |   150
        Deteriorated coating or light slime |   10–20      |  100   |   300
        Heavy slime                         |   30         |  300   |   600
        Small calcareous fouling or weed    |   40–60      |  1000  |   1000
        Medium calcareous fouling           |   70–80      |  3000  |   3000
        Heavy calcareous fouling            |   90–100     |  10000 |   10000
    \endtable


Usage
    \table
        Property            | Description             | Required    | Default value
        roughnessHeight     | Equivalent sand-grain roughness height | yes     |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type                    nutkDemirelRoughWallFunction;
        roughnessHeight         0.0004;
        value                   uniform 0.0001;
    }
    \endverbatim

Reference
    \verbatim
        Demirel, Y. K., Turan, O., & Incecik, A. (2017).
        Predicting the effect of biofouling on ship resistance
        using CFD. Applied ocean research, 62, 100-118.
    \endverbatim

See also
    Foam::nutkDemirelRoughWallFunctionFvPatchScalarField

SourceFiles
    nutkDemirelRoughWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutkDemirelRoughWallFunctionFvPatchScalarField_H
#define nutkDemirelRoughWallFunctionFvPatchScalarField_H

#include "derivedFvPatchFields/wallFunctions/nutWallFunctions/nutkWallFunction/nutkWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class nutkDemirelRoughWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutkDemirelRoughWallFunctionFvPatchScalarField
:
    public nutkWallFunctionFvPatchScalarField
{
protected:

    // Protected Member Functions

        //- Compute the roughness function
        virtual scalar E(const scalar KsPlus) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const override;


public:

    //- Runtime type information
    TypeName("nutkDemirelRoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutkDemirelRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutkDemirelRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutkDemirelRoughWallFunctionFvPatchScalarField
        //  onto a new patch
        nutkDemirelRoughWallFunctionFvPatchScalarField
        (
            const nutkDemirelRoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutkDemirelRoughWallFunctionFvPatchScalarField
        (
            const nutkDemirelRoughWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const override
        {
            return tmp<fvPatchScalarField>
            (
                new nutkDemirelRoughWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutkDemirelRoughWallFunctionFvPatchScalarField
        (
            const nutkDemirelRoughWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchScalarField>
            (
                new nutkDemirelRoughWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Return U+ based on yPlus
        scalar getUPlus(const scalar& yPlus) const override;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&) override;

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            ) override;

        // I-O

            //- Write
            virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
