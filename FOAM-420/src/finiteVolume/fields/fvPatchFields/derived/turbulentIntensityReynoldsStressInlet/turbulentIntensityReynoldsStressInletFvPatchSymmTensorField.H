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
    (c) 2017 Esi Ltd.

Class
    Foam::turbulentIntensityReynoldsStressInletFvPatchSymmTensorField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides an inlet condition for Reynolds Stress models,
    based on user-supplied turbulence intensity, defined as a fraction of the
    mean velocity:

        \f[
            R_ij = 2/3 k = 2/3*(3/2*(I |U|)^2)*delta_ij
        \f]

    where

    \vartable
        R_ij=j  | Reynolds Stress Tensor
        I       | turbulence intensity
        U       | velocity field
    \endvartable

    In the event of reverse flow, a zero-gradient condition is applied.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        intensity    | fraction of mean field [0-1] | yes    |
        U            | velocity field name     | no          | U
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type        turbulentIntensityReynoldsStressInlet;
        intensity   0.05;           // 5% turbulence
        value       uniform (1 1 1 1 1 1 );      // placeholder
    }
    \endverbatim

SeeAlso
    Foam::inletOutletFvPatchField

SourceFiles
    turbulentIntensityReynoldsStressInletFvPatchSymmTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentIntensityReynoldsStressInletFvPatchSymmTensorField_H
#define turbulentIntensityReynoldsStressInletFvPatchSymmTensorField_H

#include "fields/fvPatchFields/derived/inletOutlet/inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
 Class turbulentIntensityReynoldsStressInletFvPatchSymmTensorField Declaration
\*---------------------------------------------------------------------------*/

class turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
:
    public inletOutletFvPatchSymmTensorField
{
    // Private data

        //- Turbulent intensity as fraction of mean velocity
        scalar intensity_;

        //- Name of the velocity field
        word UName_;

        //- Name of the density field
        word rhoName_;

public:

    //- Runtime type information
    TypeName("turbulentIntensityReynoldsStressInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
        (
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
        (
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
        //  onto a new patch
        turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
        (
            const turbulentIntensityReynoldsStressInletFvPatchSymmTensorField&,
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
        (
            const turbulentIntensityReynoldsStressInletFvPatchSymmTensorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchSymmTensorField> clone() const
        {
            return tmp<fvPatchSymmTensorField>
            (
                new turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
        (
            const turbulentIntensityReynoldsStressInletFvPatchSymmTensorField&,
            const DimensionedField<symmTensor, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchSymmTensorField> clone
        (
            const DimensionedField<symmTensor, volMesh>& iF
        ) const
        {
            return tmp<fvPatchSymmTensorField>
            (
                new turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
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

        virtual void operator=(const fvPatchField<symmTensor>& pvf);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
