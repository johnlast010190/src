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
    (c) 2016 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/fixedFluxExtrapolatedPressure/fixedFluxExtrapolatedPressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "finiteVolume/fvc/fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedFluxExtrapolatedPressureFvPatchScalarField::
fixedFluxExtrapolatedPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedFluxPressureFvPatchScalarField(p, iF)
{}


Foam::fixedFluxExtrapolatedPressureFvPatchScalarField::
fixedFluxExtrapolatedPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedFluxPressureFvPatchScalarField(p, iF, dict)
{}


Foam::fixedFluxExtrapolatedPressureFvPatchScalarField::
fixedFluxExtrapolatedPressureFvPatchScalarField
(
    const fixedFluxExtrapolatedPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedFluxPressureFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::fixedFluxExtrapolatedPressureFvPatchScalarField::
fixedFluxExtrapolatedPressureFvPatchScalarField
(
    const fixedFluxExtrapolatedPressureFvPatchScalarField& wbppsf
)
:
    fixedFluxPressureFvPatchScalarField(wbppsf)
{}


Foam::fixedFluxExtrapolatedPressureFvPatchScalarField::
fixedFluxExtrapolatedPressureFvPatchScalarField
(
    const fixedFluxExtrapolatedPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedFluxPressureFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedFluxExtrapolatedPressureFvPatchScalarField::evaluate
(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    //this should be cached to avoid big overheads
    //the extrapolatedGradient switch allows the gradient to be calculated
    //as if the gradient is the same on the boundary as in the cell centre
    const tmp<volVectorField> gradp
    (
        fvc::grad
        (
            db().lookupObject<volScalarField>(internalField().name())
        )
    );

    const tmp<vectorField> dPp
    (
        gradp->boundaryField()[patch().index()].patchInternalField()
    );

    //this does not match skew correction used for grad(p)
    //how do we know if skew correction is being used or not when setting
    //the value? For now we use delta.

    //const tmp<vectorField> delta(patch().Cf() - patch().Cn());

    Field<scalar>::operator=
    (
        this->patchInternalField() + (dPp & patch().delta())
    );

    fvPatchField<scalar>::evaluate();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedFluxExtrapolatedPressureFvPatchScalarField
    );
}


// ************************************************************************* //
