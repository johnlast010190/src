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

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/translatingWallVelocity/translatingWallVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    U_(nullptr)
{}


Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    U_(Function1<vector>::New("U", dict))
{
    // Evaluate the wall velocity
    updateCoeffs();
}


Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const translatingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    U_(ptf.U_, false)
{}


Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const translatingWallVelocityFvPatchVectorField& twvpvf
)
:
    fixedValueFvPatchField<vector>(twvpvf),
    U_(twvpvf.U_, false)
{}


Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const translatingWallVelocityFvPatchVectorField& twvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(twvpvf, iF),
    U_(twvpvf.U_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::translatingWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    const vector U = U_->value(t);

    // Remove the component of U normal to the wall in case the wall is not flat
    const vectorField n(patch().nf());
    vectorField::operator=(U - n*(n & U));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::translatingWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    U_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        translatingWallVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
