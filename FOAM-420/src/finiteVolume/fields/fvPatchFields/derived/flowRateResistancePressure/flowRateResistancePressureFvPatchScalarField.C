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

#include "fields/fvPatchFields/derived/flowRateResistancePressure/flowRateResistancePressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateResistancePressureFvPatchScalarField::
flowRateResistancePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    p0_(0.0),
    R_(0.0),
    alpha_(1.0)
{
}


Foam::flowRateResistancePressureFvPatchScalarField::
flowRateResistancePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    p0_(dict.lookupOrDefault<scalar>("p0", 0.0)),
    R_(readScalar(dict.lookup("R"))),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 1.0))
{
    // Set to p0 if value is not available
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(p0_);
    }
}


Foam::flowRateResistancePressureFvPatchScalarField::
flowRateResistancePressureFvPatchScalarField
(
    const flowRateResistancePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    p0_(ptf.p0_),
    R_(ptf.R_),
    alpha_(ptf.alpha_)
{}


Foam::flowRateResistancePressureFvPatchScalarField::
flowRateResistancePressureFvPatchScalarField
(
    const flowRateResistancePressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    p0_(ptf.p0_),
    R_(ptf.R_),
    alpha_(ptf.alpha_)
{}


Foam::flowRateResistancePressureFvPatchScalarField::
flowRateResistancePressureFvPatchScalarField
(
    const flowRateResistancePressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    p0_(ptf.p0_),
    R_(ptf.R_),
    alpha_(ptf.alpha_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar Foam::flowRateResistancePressureFvPatchScalarField::pLumpedModel
(
    scalar Q
)
{
    return p0_ + R_*Q;
}

void Foam::flowRateResistancePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchVectorField& Up =
            patch().lookupPatchFieldInDb<volVectorField, vector>(db(), "U");

    const vectorField& Sf = patch().Sf();

    // Patch normal velocity
    // Obs.: the patch field is used to take advantage of relaxation factors of
    // the SIMPLE algorithm.
    scalarField Upn( (Up & Sf)/patch().magSf() );

    // Calculate current flow rate
    scalar Q = gSum(Upn*patch().magSf());

    // Total patch area
    scalar patchArea = gSum(patch().magSf());

    // Calculate dynamic head for inflow
    scalarField Kin( alpha_*neg(Upn)*0.5*mag(Upn)*Upn );
    scalar Kmean = gSum(Kin*patch().magSf())/patchArea;

    // Calculate equivalent pressure to current flow rate
    scalar p = pLumpedModel(Q);

    if (debug)
    {
        Info<< "* " << patch().name() << ":  "
             << "Q = " << Q << ", "
             << "p = " << p << ", "
             << "Kmean = " << Kmean << endl;
    }

    fixedValueFvPatchScalarField::forceAssign(p - Kmean + Kin);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::flowRateResistancePressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<scalar>(os, "p0", 0.0, p0_);
    os.writeEntry("R", R_);
    writeEntryIfDifferent<scalar>(os, "alpha", 1.0, alpha_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        flowRateResistancePressureFvPatchScalarField
    );
}

// ************************************************************************* //
