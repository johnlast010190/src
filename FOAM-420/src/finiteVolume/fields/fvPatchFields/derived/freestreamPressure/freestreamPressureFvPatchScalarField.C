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
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/freestreamPressure/freestreamPressureFvPatchScalarField.H"
#include "fields/fvPatchFields/derived/freestreamVelocity/freestreamVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
//
Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    UName_("U"),
    supersonic_(false)
{
}

Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    supersonic_
    (
        dict.lookupOrDefault<Switch>("supersonic", false)
    )
{
    referenceFrameFvPatch<scalar>::read(dict);
    this->patchType() = dict.lookupOrDefault<word>("patchType", word::null);
    freestreamValue() = scalarField("freestreamValue", dict, p.size());

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(freestreamValue());
    }

    refGrad() = Zero;
    if (dict.found("valueFraction"))
    {
        valueFraction()  = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        valueFraction() = 0;
    }
}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    UName_(psf.UName_),
    supersonic_(psf.supersonic_)
{}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& psf
)
:
    mixedFvPatchScalarField(psf),
    UName_(psf.UName_),
    supersonic_(psf.supersonic_)
{}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    UName_(psf.UName_),
    supersonic_(psf.supersonic_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freestreamPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    if (!isA<freestreamVelocityFvPatchVectorField>(Up))
    {
        FatalErrorInFunction
        << "Cannot use freestreamPressure without freestreamVelocity for U"
        << nl << abort(FatalError);
    }
    const freestreamVelocityFvPatchVectorField& Ustp =
        refCast<const freestreamVelocityFvPatchVectorField>(Up);

    Field<scalar>& vf = valueFraction();
    if (supersonic_)
    {
        vf = Ustp.valueFraction();
    }
    else
    {
        vf = 1-Ustp.valueFraction();
    }

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::freestreamPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    freestreamValue().writeEntry("freestreamValue", os);
    os.writeEntry("supersonic", supersonic_);
    valueFraction().writeEntry("valueFraction", os);
    referenceFrameFvPatch<scalar>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        freestreamPressureFvPatchScalarField
    );
}

// ************************************************************************* //
