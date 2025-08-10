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
    (c) 2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "pressureVibratingWindowFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pressureVibratingWindowFvPatchScalarField::pressureVibratingWindowFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    accelerationName_ (word("w")),
    rho_            (scalar(1.205))
{}


Foam::pressureVibratingWindowFvPatchScalarField::pressureVibratingWindowFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    accelerationName_ (dict.lookup("accelerationName")),
    rho_            (readScalar(dict.lookup("rho")))

{
    this->forceAssign(this->patchInternalField()
	  + this->gradient()/this->patch().deltaCoeffs());
}


Foam::pressureVibratingWindowFvPatchScalarField::pressureVibratingWindowFvPatchScalarField
(
    const pressureVibratingWindowFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    accelerationName_ (ptf.accelerationName_),
    rho_            (ptf.rho_)
{}


Foam::pressureVibratingWindowFvPatchScalarField::pressureVibratingWindowFvPatchScalarField
(
    const pressureVibratingWindowFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    accelerationName_ (ptf.accelerationName_),
    rho_            (ptf.rho_)
{}


Foam::pressureVibratingWindowFvPatchScalarField::pressureVibratingWindowFvPatchScalarField
(
    const pressureVibratingWindowFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    accelerationName_ (ptf.accelerationName_),
    rho_            (ptf.rho_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::pressureVibratingWindowFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const areaScalarField& ac = this->db().objectRegistry::lookupObject<areaScalarField>(accelerationName_);
    this->gradient() = - (rho_ * ac.primitiveField());

    // Sets Updated to true
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::pressureVibratingWindowFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);

    os.writeEntry("accelerationName", accelerationName_);
    os.writeEntry("rho", rho_);

    this->writeEntry("value", os);
}


// ************************************************************************* //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        pressureVibratingWindowFvPatchScalarField
    );
}
// ************************************************************************* //
