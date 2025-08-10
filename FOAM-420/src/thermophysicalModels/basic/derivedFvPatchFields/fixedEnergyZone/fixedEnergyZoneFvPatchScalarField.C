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

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "basicThermo/basicThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "derivedFvPatchFields/fixedEnergyZone/fixedEnergyZoneFvPatchScalarField.H"
#include "fields/fvPatchFields/derived/fixedValueZone/fixedValueZoneFvPatchField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedEnergyZoneFvPatchScalarField::
fixedEnergyZoneFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedEnergyFvPatchScalarField(p, iF)
{}


Foam::fixedEnergyZoneFvPatchScalarField::
fixedEnergyZoneFvPatchScalarField
(
    const fixedEnergyZoneFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedEnergyFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::fixedEnergyZoneFvPatchScalarField::
fixedEnergyZoneFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedEnergyFvPatchScalarField(p, iF, dict)
{}


Foam::fixedEnergyZoneFvPatchScalarField::
fixedEnergyZoneFvPatchScalarField
(
    const fixedEnergyZoneFvPatchScalarField& tppsf
)
:
    fixedEnergyFvPatchScalarField(tppsf)
{}


Foam::fixedEnergyZoneFvPatchScalarField::
fixedEnergyZoneFvPatchScalarField
(
    const fixedEnergyZoneFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedEnergyFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedEnergyZoneFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fixedEnergyFvPatchScalarField::updateCoeffs();
}


void Foam::fixedEnergyZoneFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    const polyPatch& pthis = this->patch().patch();
    const polyMesh& pMesh =  pthis.boundaryMesh().mesh();
    const label patchi = patch().index();

    const indirectPolyPatch& gibPolyPatch =
        refCast<const indirectPolyPatch>(pthis);

    const word faceZoneName = pMesh.faceZones()[gibPolyPatch.zoneId()].name();

    const label& zoneId = pMesh.cellZones().findZoneID
    (
        "inactive_"+faceZoneName
    );

    const basicThermo& thermo = basicThermo::lookupThermo(*this);

    const fixedValueZoneFvPatchField<scalar>& Tw =
        refCast<const fixedValueZoneFvPatchField<scalar>>
        (
            thermo.T().boundaryField()[patchi]
        );

    const fixedValueZoneFvPatchField<scalar>& pw =
        refCast<const fixedValueZoneFvPatchField<scalar>>
        (
            thermo.p().boundaryField()[patchi]
        );

    const labelList& zoneCells = pMesh.cellZones()[zoneId];
    scalarField pField(zoneCells.size(), pw.cellZoneValue());
    scalarField TField(zoneCells.size(), Tw.cellZoneValue());

    tmp<Foam::scalarField> tE = thermo.he
    (
        pField,
        TField,
        zoneCells
    );

    matrix.setValues
    (
        zoneCells,
        tE()
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedEnergyZoneFvPatchScalarField
    );
}

// ************************************************************************* //
