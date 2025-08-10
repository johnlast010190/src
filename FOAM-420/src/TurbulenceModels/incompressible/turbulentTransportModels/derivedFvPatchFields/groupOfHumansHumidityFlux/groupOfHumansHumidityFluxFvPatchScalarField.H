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
    (c) 2011-2013 OpenFOAM Foundation

Class
    Foam::groupOfHumansHumidityFluxFvPatchScalarField

Description
    Empirical mass flux boundary for water vapour.

SourceFiles
    groupOfHumansHumidityFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef groupOfHumansHumidityFluxFvPatchScalarFields_H
#define groupOfHumansHumidityFluxFvPatchScalarFields_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "derivedFvPatchFields/turbulentMassFluxConcentration/turbulentMassFluxConcentrationFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
     Class groupOfHumansHumidityFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class groupOfHumansHumidityFluxFvPatchScalarField
:
    public turbulentMassFluxConcentrationFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("incompressible::groupOfHumansHumidityFlux");


    // Constructors

        //- Construct from patch and internal field
        groupOfHumansHumidityFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        groupOfHumansHumidityFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  groupOfHumansHumidityFluxFvPatchScalarField onto
        //  a new patch
        groupOfHumansHumidityFluxFvPatchScalarField
        (
            const groupOfHumansHumidityFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        groupOfHumansHumidityFluxFvPatchScalarField
        (
            const groupOfHumansHumidityFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new groupOfHumansHumidityFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        groupOfHumansHumidityFluxFvPatchScalarField
        (
            const groupOfHumansHumidityFluxFvPatchScalarField&,
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
                new groupOfHumansHumidityFluxFvPatchScalarField
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
