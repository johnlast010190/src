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
    (c) 2018 Esi Ltd.

Class
    Foam::fixedEnergyZoneFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition is used in the slave part of the GIBs to deactivate
    the solution.
    It is constructed on the fly.
    Check basicThermo.C

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedEnergyZone;
        value           uniform 0;
        cellZoneValue   0;
    }
    \endverbatim

See also
    Foam::fixedEnergyZoneFvPatchScalarField

SourceFiles
    fixedEnergyZoneFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedEnergyZoneFvPatchScalarField_H
#define fixedEnergyZoneFvPatchScalarField_H

#include "derivedFvPatchFields/fixedEnergy/fixedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class fixedEnergyZoneFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedEnergyZoneFvPatchScalarField
:
    public fixedEnergyFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("fixedEnergyZone");


    // Constructors

        //- Construct from patch and internal field
        fixedEnergyZoneFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedEnergyZoneFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedEnergyZoneFvPatchScalarField
        // onto a new patch
        fixedEnergyZoneFvPatchScalarField
        (
            const fixedEnergyZoneFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedEnergyZoneFvPatchScalarField
        (
            const fixedEnergyZoneFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedEnergyZoneFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedEnergyZoneFvPatchScalarField
        (
            const fixedEnergyZoneFvPatchScalarField&,
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
                new fixedEnergyZoneFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Manipulate matrix with given weights
            virtual void manipulateMatrix
            (
                fvMatrix<scalar>&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
