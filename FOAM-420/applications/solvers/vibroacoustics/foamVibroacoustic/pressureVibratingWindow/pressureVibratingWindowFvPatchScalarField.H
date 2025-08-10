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
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::pressureVibratingWindowFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition supplies a fixed gradient condition based on the
    window accelaration, such that the patch values are calculated using:

        \f[
            x_p = x_c - \rho \frac{\partial w}{\partial t}
        \f]

    where
    \vartable
        x_p      | patch values
        x_c      | internal field values
        \rho     | air density
        a        | window acceleration
    \endvartable

Usage
    \table
        Property        | Description              | Required    | Default value
        accelerationName  | name of acceleration field | yes         |
        rho             | density value            | yes         |
        gradient        | gradient                 | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            pressureVibratingWindow;
        accelerationName  a;
        rho             1.205;
        gradient        uniform 0;
    }
    \endverbatim

SourceFiles
    pressureVibratingWindowFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureVibratingWindowFvPatchScalarField_H
#define pressureVibratingWindowFvPatchScalarField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/fixedGradient/fixedGradientFvPatchFields.H"
#include "include/faCFD.H"
#include "finiteArea/fac/facD2dt2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class pressureVibratingWindowFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class pressureVibratingWindowFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    word    accelerationName_;
    scalar  rho_;

public:

    //- Runtime type information
    TypeName("pressureVibratingWindow");


    // Constructors

        //- Construct from patch and internal field
        pressureVibratingWindowFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureVibratingWindowFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given pressureVibratingWindowFvPatchScalarField
        //  onto a new patch
        pressureVibratingWindowFvPatchScalarField
        (
            const pressureVibratingWindowFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureVibratingWindowFvPatchScalarField
        (
            const pressureVibratingWindowFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new pressureVibratingWindowFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureVibratingWindowFvPatchScalarField
        (
            const pressureVibratingWindowFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new pressureVibratingWindowFvPatchScalarField(*this, iF)
            );
        }


    // Member functions
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
