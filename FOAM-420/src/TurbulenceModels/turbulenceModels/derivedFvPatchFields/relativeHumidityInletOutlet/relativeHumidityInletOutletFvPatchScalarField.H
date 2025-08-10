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
    (c) 2010-2011, Esi Ltd

Class
    Foam::relativeHumidityInletOutletFvPatchScalarField

Description
    Calculate humidity ratio from specified relative humidity

    Example of the boundary condition specification:
    @verbatim
        inlet
        {
            type        relativeHumidityInletOutlet;
            relativeHumidity  uniform 0.5;    // 50% humidity
            T           T;
            p           p;
            value       uniform 0.001;      // placeholder
        rho         rho;  //rhoEff;
        }
    @endverbatim

SourceFiles
    relativeHumidityInletOutletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef relativeHumidityInletOutletFvPatchScalarField_H
#define relativeHumidityInletOutletFvPatchScalarField_H

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
          Class relativeHumidityInletOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class relativeHumidityInletOutletFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{
    // Private data

        //- Turbulent intensity as fraction of mean velocity
        scalarField relativeHumidity_;

        word phiName_;
        word TName_;
        word PName_;
        word rhoName_;

        // physical constants

            //- Molecular weight of vapour phase
            scalar Mvap_;

            //- Molecular weight of mixture
            scalar Mmix_;

private:

    //- lookup patch density
    tmp<scalarField> patchDensity();

public:

    //- Runtime type information
    TypeName("relativeHumidityInletOutlet");


    // Constructors

        //- Construct from patch and internal field
        relativeHumidityInletOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        relativeHumidityInletOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  relativeHumidityInletOutletFvPatchScalarField
        //  onto a new patch
        relativeHumidityInletOutletFvPatchScalarField
        (
            const relativeHumidityInletOutletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        relativeHumidityInletOutletFvPatchScalarField
        (
            const relativeHumidityInletOutletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new relativeHumidityInletOutletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        relativeHumidityInletOutletFvPatchScalarField
        (
            const relativeHumidityInletOutletFvPatchScalarField&,
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
                new relativeHumidityInletOutletFvPatchScalarField
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


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

    // Member operators

        virtual void operator=(const fvPatchField<scalar>& pvf);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
