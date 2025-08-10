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
    Foam::nutUTabulatedWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions.  As input, the user specifies a look-up table
    of U+ as a function of near-wall Reynolds number.

Usage
    \table
        Property     | Description             | Required    | Default value
        uPlusTable   | U+ as a function of Re file name | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutTabulatedWallFunction;
        uPlusTable      myUPlusTable;
    }
    \endverbatim

Note
    The tables are not registered since the same table object may be used for
    more than one patch.

See also
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    nutUTabulatedWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutUTabulatedWallFunctionFvPatchScalarField_H
#define nutUTabulatedWallFunctionFvPatchScalarField_H

#include "derivedFvPatchFields/wallFunctions/nutWallFunctions/nutWallFunction/nutWallFunctionFvPatchScalarField.H"
#include "interpolations/interpolationTable/interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class nutUTabulatedWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUTabulatedWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{

protected:

    // Protected data

        //- FileName
        word fileName_;

        //- File containing u+(y+) data
        fileName uPlusTableFileData_;

        //- U+ table
        interpolationTable<scalar> uPlusTable_;


    // Protected Member Functions

        fileName getTableFileName(const dictionary& dict);

        //- Calculate the turbulence viscosity
        tmp<scalarField> calcNut() const override;

        //- Calculate wall u+ from table
        virtual tmp<scalarField> calcuTau
        (
            const scalarField& yPlus,
            const scalarField& Uw
        ) const;


public:

    //- Runtime type information
    TypeName("nutUTabulatedWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutUTabulatedWallFunctionFvPatchScalarField
        //  onto a new patch
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const nutUTabulatedWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const nutUTabulatedWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        tmp<fvPatchScalarField> clone() const override
        {
            return tmp<fvPatchScalarField>
            (
                new nutUTabulatedWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const nutUTabulatedWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchScalarField>
            (
                new nutUTabulatedWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            tmp<scalarField> yPlus() const override;

        // I-O

            //- Write
             void write(Ostream& os) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
