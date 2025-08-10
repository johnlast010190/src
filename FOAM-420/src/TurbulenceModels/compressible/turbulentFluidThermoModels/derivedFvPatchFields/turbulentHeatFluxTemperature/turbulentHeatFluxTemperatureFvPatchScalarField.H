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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2022 Esi Ltd

Class
    Foam::turbulentHeatFluxTemperatureFvPatchScalarField

Description
    Fixed heat boundary condition to specify temperature gradient. Input
    heat source either specified in terms of an absolute power [W], or as a
    flux [W/m^2].

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        heatSource   | 'power' [W] or 'flux' [W/m^2] | yes |
        q            | heat power or flux field      | yes |
        kappaMethod  | inherited from boundaryKappa | yes |
        kappa         | inherited from boundaryKappa | yes |
        kappaName    | inherited from boundaryKappa | yes |
        Qr           | name of the radiative flux field | yes |
        value        | initial temperature value | no | calculated
        gradient     | initial gradient value | no | 0.0
    qtime        | heat power or flux field in tabular format
    \endtable

    Note: If needed, both 'value' and 'gradient' must be defined to be used.

    Example usage:
    \verbatim
    hotWall
    {
        type            compressible::turbulentHeatFluxTemperature;
        heatSource      flux;
        q               uniform 10;
        kappaMethod     fluidThermo;
        kappa            none;
        Qr              none;
    }
    \endverbatim

    For time dependent inputs

    \verbatim
    hotWall
    {
        type            compressible::turbulentHeatFluxTemperature;
    heatSource      fluxTV;
        kappa           fluidThermo;
    kappaName       none;
        outOfBounds     clamp;
        qtime           table ( ( 0.0 0.0 ) ( 10.0 100.0 ) ( 20.0 200.0 ) );
    }
    \endverbatim


SeeAlso
    Foam::boundaryKappa

SourceFiles
    turbulentHeatFluxTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentHeatFluxTemperatureFvPatchScalarFields_H
#define turbulentHeatFluxTemperatureFvPatchScalarFields_H

#include "fields/fvPatchFields/basic/fixedGradient/fixedGradientFvPatchFields.H"
#include "turbulentFluidThermoModels/derivedFvPatchFields/boundaryKappa/boundaryKappa.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
      Class turbulentHeatFluxTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentHeatFluxTemperatureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField,
    public boundaryKappa
{
public:

    // Data types

        //- Enumeration listing the possible hest source input modes
        enum heatSourceType
        {
            hsPower,
            hsFlux,
        hsTVPower,
        hsTVFlux
        };


private:

    // Private data

        //- Heat source type names
        static const NamedEnum<heatSourceType, 4> heatSourceTypeNames_;

        //- Heat source type
        heatSourceType heatSource_;


protected:

        //- Heat power [W] or flux [W/m2]
        scalarField q_;

        //- Time Varing Flux
        autoPtr<Function1<scalar>> qtimevar_;

        //- Name of radiative in flux field
        word QrName_;


public:

    //- Runtime type information
    TypeName("compressible::turbulentHeatFluxTemperature");


    // Constructors

        //- Construct from patch and internal field
        turbulentHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentHeatFluxTemperatureFvPatchScalarField onto
        //  a new patch
        turbulentHeatFluxTemperatureFvPatchScalarField
        (
            const turbulentHeatFluxTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentHeatFluxTemperatureFvPatchScalarField
        (
            const turbulentHeatFluxTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentHeatFluxTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentHeatFluxTemperatureFvPatchScalarField
        (
            const turbulentHeatFluxTemperatureFvPatchScalarField&,
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
                new turbulentHeatFluxTemperatureFvPatchScalarField
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

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
