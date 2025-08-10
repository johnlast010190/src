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
    (c) 2022 Esi Ltd.

Class
    Foam::mappedTotalTemperatureProfileFvPatchScalarField

Description
    This boundary condition provides a total temperature condition. It maps
    total temperature values from a file to a patch (using the base class), and
    then computes static temperature values.

Usage
    \table
        Property     | Description             | Required    | Default value
        U            | Velocity field name     | no          | U
        phi          | Flux field name         | no          | phi
        psi          | Compressibility field name | no       | thermo:psi
        gamma        | ratio of specific heats (Cp/Cv) | yes |
        T0           | total temperature       | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mappedTotalTemperatureProfile;
        file            "data.csv";
        gamma           1.4;
        T0              uniform 800;
        field           Ttot;
    }
    \endverbatim


See also
    Foam::mappedProfileFvPatchField<Type>

SourceFiles
    mappedTotalTemperatureProfileFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedTotalTemperatureProfileFvPatchScalarField_H
#define mappedTotalTemperatureProfileFvPatchScalarField_H

#include "fields/fvPatchFields/derived/mappedProfile/mappedProfileFvPatchFields.H"
#include "fields/fvPatchFields/derived/mappedProfile/mappedProfileFvPatchField.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class mappedTotalTemperatureProfileFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mappedTotalTemperatureProfileFvPatchScalarField
:
    public mappedProfileFvPatchField<scalar>
{

    // Private data

        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the compressibility field used to calculate the wave speed
        word psiName_;

        //- Heat capacity ratio
        scalar gamma_;

    // Private functions

        virtual void updateValues();


public:

    //- Runtime type information
    TypeName("mappedTotalTemperatureProfile");


    // Constructors

        //- Construct from patch and internal field
        mappedTotalTemperatureProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedTotalTemperatureProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mappedTotalTemperatureProfileFvPatchScalarField
        //  onto a new patch
        mappedTotalTemperatureProfileFvPatchScalarField
        (
            const mappedTotalTemperatureProfileFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mappedTotalTemperatureProfileFvPatchScalarField
        (
            const mappedTotalTemperatureProfileFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mappedTotalTemperatureProfileFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mappedTotalTemperatureProfileFvPatchScalarField
        (
            const mappedTotalTemperatureProfileFvPatchScalarField&,
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
                new mappedTotalTemperatureProfileFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
