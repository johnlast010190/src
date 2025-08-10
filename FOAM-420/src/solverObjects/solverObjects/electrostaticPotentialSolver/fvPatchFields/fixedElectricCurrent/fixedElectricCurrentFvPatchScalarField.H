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
    (c) 2019-2023 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::fixedElectricCurrentFvPatchScalarField

Group
    grpThermoBoundaryConditions grpWallBoundaryConditions

Description
    This boundary condition applies a prescribed electric current
    on an external wall.

    If no magnetic field is present:

    \f[
        J = \sigma E
    \f]

    where \f$ J \f$ is the current density and \f$ E \f$ the electric field.

    The electric field is given by

    \f[
        E = \grad V
    \f]


    The thermal conductivity \c sigma is looked up from the database.

Usage
    \table
    Property     | Description                 | Required | Default value
    I            | Electric current [A] (supports time-varying Function1)    | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedElectricCurrent;
        I               2.0;
        value           $internalField;
    }
    \endverbatim

See also
    Foam::fv::electrostaticPotentialSolver
    Foam::fixedGradientFvPatchScalarField

SourceFiles
    fixedElectricCurrentFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedElectricCurrentFvPatchScalarField_H
#define fixedElectricCurrentFvPatchScalarField_H

#include "fields/fvPatchFields/basic/fixedGradient/fixedGradientFvPatchFields.H"
#include "fvPatchFields/electricalBoundaryBase/electricalBoundaryBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class fixedElectricCurrentFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedElectricCurrentFvPatchScalarField
:
    public fixedGradientFvPatchScalarField,
    public electricalBoundaryBase
{

private:

    // Private data

        //- Electric current [A]
        autoPtr<Function1<scalar>> I_;


public:

    //- Runtime type information
    TypeName("fixedElectricCurrent");


    // Constructors

        //- Construct from patch and internal field
        fixedElectricCurrentFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedElectricCurrentFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedElectricCurrentFvPatchScalarField
        //  onto a new patch
        fixedElectricCurrentFvPatchScalarField
        (
            const fixedElectricCurrentFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedElectricCurrentFvPatchScalarField
        (
            const fixedElectricCurrentFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedElectricCurrentFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedElectricCurrentFvPatchScalarField
        (
            const fixedElectricCurrentFvPatchScalarField&,
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
                new fixedElectricCurrentFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
