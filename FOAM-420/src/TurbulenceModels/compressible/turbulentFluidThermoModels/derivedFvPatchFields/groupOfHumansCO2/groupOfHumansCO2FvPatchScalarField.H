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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2009 OpenCFD Ltd.

Class
    Foam::groupOfHumansCO2FvPatchScalarField

Description
    Empirical mass flux boundary for CO2.

SourceFiles
    groupOfHumansCO2FvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef groupOfHumansCO2FvPatchScalarFields_H
#define groupOfHumansCO2FvPatchScalarFields_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "derivedFvPatchFields/turbulentMassFluxConcentration/turbulentMassFluxConcentrationFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
     Class groupOfHumansCO2FvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class groupOfHumansCO2FvPatchScalarField
:
    public turbulentMassFluxConcentrationFvPatchScalarField
{

private:

    //- input data

        // CO2 gas constant
        scalar RCO2_;

public:

    //- Runtime type information
    TypeName("compressible::groupOfHumansCO2");


    // Constructors

        //- Construct from patch and internal field
        groupOfHumansCO2FvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        groupOfHumansCO2FvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  groupOfHumansCO2FvPatchScalarField onto
        //  a new patch
        groupOfHumansCO2FvPatchScalarField
        (
            const groupOfHumansCO2FvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        groupOfHumansCO2FvPatchScalarField
        (
            const groupOfHumansCO2FvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new groupOfHumansCO2FvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        groupOfHumansCO2FvPatchScalarField
        (
            const groupOfHumansCO2FvPatchScalarField&,
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
                new groupOfHumansCO2FvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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
