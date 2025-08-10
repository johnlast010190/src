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
    (c) 2006-2009 OpenCFD Ltd.
    (c) 2010-2012 Esi Ltd.

Class
    Foam::nutTurbulentIntensityLengthScaleInletFvPatchScalarField

Description
    Calculate turbulent viscosity from the provided intensity and length scale.
    Useful for nuTilda on inlets

    Example of the boundary condition specification:
    @verbatim
        inlet
        {
            type        nutTurbulentIntensityLengthScaleInlet;
            intensity   0.05;           // 5% turbulence
            length      0.01;           // 1 cm
            value       uniform 1;      // placeholder
        }
    @endverbatim

SourceFiles
    nutTurbulentIntensityLengthScaleInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutTurbulentIntensityLengthScaleInletFvPatchScalarField_H
#define nutTurbulentIntensityLengthScaleInletFvPatchScalarField_H

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutTurbulentIntensityLengthScaleInletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nutTurbulentIntensityLengthScaleInletFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{
    // Private data

        //- Turbulent intensity as fraction of mean velocity
        scalar intensity_;

        //- Turbulent length scale
        scalar length_;

        //- Cmu
        scalar Cmu_;

        word phiName_;

        //- U field name
        word UName_;

public:

    //- Runtime type information
    TypeName("nutTurbulentIntensityLengthScaleInlet");


    // Constructors

        //- Construct from patch and internal field
        nutTurbulentIntensityLengthScaleInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutTurbulentIntensityLengthScaleInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutTurbulentIntensityLengthScaleInletFvPatchScalarField
        //  onto a new patch
        nutTurbulentIntensityLengthScaleInletFvPatchScalarField
        (
            const nutTurbulentIntensityLengthScaleInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutTurbulentIntensityLengthScaleInletFvPatchScalarField
        (
            const nutTurbulentIntensityLengthScaleInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutTurbulentIntensityLengthScaleInletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        nutTurbulentIntensityLengthScaleInletFvPatchScalarField
        (
            const nutTurbulentIntensityLengthScaleInletFvPatchScalarField&,
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
                new nutTurbulentIntensityLengthScaleInletFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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
