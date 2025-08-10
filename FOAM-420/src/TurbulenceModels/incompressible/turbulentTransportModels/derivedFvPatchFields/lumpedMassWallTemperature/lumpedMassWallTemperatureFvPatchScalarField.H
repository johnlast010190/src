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
    (c) 2020 Esi Ltd.
    (c) 2016 OpenCFD Ltd.

Class
    lumpedMassWallTemperatureFvPatchScalarField

Description
    Employs a lumped mass model for temperature.

    It considers a single temperature value for the whole patch and evaluates
    the temperature evolution using the net heat flux into the patch.

Usage
    \table
        Property     | Description                 | Required | Default value
        Cp           | Heat capacity [J/kg.K]      | yes      |
        mass         | Total mass [kg]             | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            incompressible::lumpedMassWallTemperature;
        mass            1000;
        Cp              4100;
        value           uniform 350.0;
    }
    \endverbatim

SourceFiles
    lumpedMassWallTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedMassWallTemperatureFvPatchScalarField_H
#define lumpedMassWallTemperatureFvPatchScalarField_H

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
         Class lumpedMassWallTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class lumpedMassWallTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
private:

    // Private data

        //- Specifc heat capacity [J/kg/K]
        scalar Cp_;

        //- Mass [kg]
        scalar mass_;

        //- Current time index (used for updating)
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("incompressible::lumpedMassWallTemperature");


    // Constructors

        //- Construct from patch and internal field
        lumpedMassWallTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        lumpedMassWallTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // lumpedMassWallTemperatureFvPatchScalarField
        // onto a new patch
        lumpedMassWallTemperatureFvPatchScalarField
        (
            const lumpedMassWallTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        lumpedMassWallTemperatureFvPatchScalarField
        (
            const lumpedMassWallTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new lumpedMassWallTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        lumpedMassWallTemperatureFvPatchScalarField
        (
            const lumpedMassWallTemperatureFvPatchScalarField&,
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
                new lumpedMassWallTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();


        // I-O

        //- Write
        void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
