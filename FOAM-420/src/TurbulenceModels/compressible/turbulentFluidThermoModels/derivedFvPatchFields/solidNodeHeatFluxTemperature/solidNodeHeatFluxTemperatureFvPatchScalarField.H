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
    (c) 2010-2019 Esi Ltd.

Class
    Foam::solidNodeHeatFluxTemperatureFvPatchScalarField

Description
    Steady state (lumped style) thermal boundary condition. The thermal
    power q specified in input is injected in a solid node. A thermal
    circuit with resistances connecting the solid node with each face of
    the patch is solved. Solid layers can be specified between node and wall.
    This boundary condition can be used to prevent the possibility of hot
    spots on patches, as a non uniform heat flux is calculated and set
    on each face of the pach.

    heatedPatch
    {
        type            compressible::solidNodeHeatFluxTemperature;
        gradient        uniform 0;
        q               constant 1500;
        value           uniform 300;
        layers
        {
            layer1
            {
                thickness       0.1;
                kappa           10;
                rcontact        0;
            }
        }
    }

SourceFiles
    solidNodeHeatFluxTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef solidNodeHeatFluxTemperatureFvPatchScalarFields_H
#define solidNodeHeatFluxTemperatureFvPatchScalarFields_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/fixedGradient/fixedGradientFvPatchFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "layersProperties/layersProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
/*---------------------------------------------------------------------------*\
     Class solidNodeHeatFluxTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class solidNodeHeatFluxTemperatureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

// Private

        //- Heat power field [W]
        autoPtr<Function1<scalar>> q_;

        // layers material properties
        layersProperties layersProp_;

public:

    //- Runtime type information
    TypeName("compressible::solidNodeHeatFluxTemperature");


    // Constructors

        //- Construct from patch and internal field
        solidNodeHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidNodeHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidNodeHeatFluxTemperatureFvPatchScalarField onto
        //  a new patch
        solidNodeHeatFluxTemperatureFvPatchScalarField
        (
            const solidNodeHeatFluxTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidNodeHeatFluxTemperatureFvPatchScalarField
        (
            const solidNodeHeatFluxTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new solidNodeHeatFluxTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidNodeHeatFluxTemperatureFvPatchScalarField
        (
            const solidNodeHeatFluxTemperatureFvPatchScalarField&,
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
                new solidNodeHeatFluxTemperatureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        // Mapping functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(const fvPatchFieldMapper&);

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchScalarField&,
            const labelList&
        );

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
