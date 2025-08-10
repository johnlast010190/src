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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/interstitialInletVelocity/interstitialInletVelocityFvPatchVectorField.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interstitialInletVelocityFvPatchVectorField::
interstitialInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    inletVelocity_(p.size(), Zero),
    alphaName_("alpha")
{}


Foam::interstitialInletVelocityFvPatchVectorField::
interstitialInletVelocityFvPatchVectorField
(
    const interstitialInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    inletVelocity_(mapper(ptf.inletVelocity_)),
    alphaName_(ptf.alphaName_)
{}


Foam::interstitialInletVelocityFvPatchVectorField::
interstitialInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    inletVelocity_("inletVelocity", dict, p.size()),
    alphaName_(dict.lookupOrDefault<word>("alpha", "alpha"))
{}


Foam::interstitialInletVelocityFvPatchVectorField::
interstitialInletVelocityFvPatchVectorField
(
    const interstitialInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    inletVelocity_(ptf.inletVelocity_),
    alphaName_(ptf.alphaName_)
{}


Foam::interstitialInletVelocityFvPatchVectorField::
interstitialInletVelocityFvPatchVectorField
(
    const interstitialInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    inletVelocity_(ptf.inletVelocity_),
    alphaName_(ptf.alphaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interstitialInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    m(inletVelocity_, inletVelocity_);
}


void Foam::interstitialInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const interstitialInletVelocityFvPatchVectorField& tiptf =
        refCast<const interstitialInletVelocityFvPatchVectorField>(ptf);

    inletVelocity_.rmap(tiptf.inletVelocity_, addr);
}


void Foam::interstitialInletVelocityFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::autoMapGIB(mapper);
    mapper.map(inletVelocity_, vector::zero);
}


void Foam::interstitialInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& alphap =
        patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), alphaName_);

    forceAssign(inletVelocity_/alphap);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::interstitialInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntryIfDifferent<word>(os, "alpha", "alpha", alphaName_);
    inletVelocity_.writeEntry("inletVelocity", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       interstitialInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
