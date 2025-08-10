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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/fixedPressureCompressibleDensity/fixedPressureCompressibleDensityFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedPressureCompressibleDensityFvPatchScalarField::
fixedPressureCompressibleDensityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    pName_("p")
{}


Foam::fixedPressureCompressibleDensityFvPatchScalarField::
fixedPressureCompressibleDensityFvPatchScalarField
(
    const fixedPressureCompressibleDensityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    pName_(ptf.pName_)
{}


Foam::fixedPressureCompressibleDensityFvPatchScalarField::
fixedPressureCompressibleDensityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    pName_(dict.lookupOrDefault<word>("p", "p"))
{}


Foam::fixedPressureCompressibleDensityFvPatchScalarField::
fixedPressureCompressibleDensityFvPatchScalarField
(
    const fixedPressureCompressibleDensityFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    pName_(ptf.pName_)
{}


Foam::fixedPressureCompressibleDensityFvPatchScalarField::
fixedPressureCompressibleDensityFvPatchScalarField
(
    const fixedPressureCompressibleDensityFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    pName_(ptf.pName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedPressureCompressibleDensityFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& pp =
        patch().lookupPatchFieldInDb<volScalarField, scalar>(this->db(), pName_);

    const dictionary& thermoProps =
        db().lookupObject<IOdictionary>("thermodynamicProperties");

    const scalar rholSat =
        dimensionedScalar("rholSat", dimDensity, thermoProps).value();

    const scalar pSat =
        dimensionedScalar("pSat", dimPressure, thermoProps).value();

    const scalar psil =
        dimensionedScalar("psil", thermoProps.lookup("psil")).value();

    forceAssign(rholSat + psil*(pp - pSat));

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::fixedPressureCompressibleDensityFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedPressureCompressibleDensityFvPatchScalarField
    );
}

// ************************************************************************* //
