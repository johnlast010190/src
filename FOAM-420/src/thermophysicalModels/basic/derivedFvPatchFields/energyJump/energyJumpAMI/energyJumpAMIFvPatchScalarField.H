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
    (c) 2012-2016 OpenFOAM Foundation

Class
    Foam::energyJumpAMIFvPatchScalarField

Group
    grpThermoBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary condition provides an energy jump condition across a pair
    of coupled patches with an arbitrary mesh interface (AMI).  It is not
    applied directly, but is employed on-the-fly when converting temperature
    boundary conditions into energy.

See also
    Foam::fixedJumpAMIFvPatchField

SourceFiles
    energyJumpAMIFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef energyJumpAMIFvPatchScalarField_H
#define energyJumpAMIFvPatchScalarField_H

#include "fields/fvPatchFields/derived/fixedJumpAMI/fixedJumpAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class energyJumpAMIFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class energyJumpAMIFvPatchScalarField
:
    public fixedJumpAMIFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("energyJumpAMI");

    // Constructors

        //- Construct from patch and internal field
        energyJumpAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        energyJumpAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given energyJumpAMIFvPatchScalarField onto a
        //  new patch
        energyJumpAMIFvPatchScalarField
        (
            const energyJumpAMIFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        energyJumpAMIFvPatchScalarField
        (
            const energyJumpAMIFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new energyJumpAMIFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        energyJumpAMIFvPatchScalarField
        (
            const energyJumpAMIFvPatchScalarField&,
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
                new energyJumpAMIFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
