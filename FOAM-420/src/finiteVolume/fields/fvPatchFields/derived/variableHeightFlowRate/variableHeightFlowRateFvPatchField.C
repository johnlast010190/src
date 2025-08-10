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
    (c) 2017 OpenCFD Ltd.
    (c) 2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/variableHeightFlowRate/variableHeightFlowRateFvPatchField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::variableHeightFlowRateFvPatchScalarField
::variableHeightFlowRateFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    lowerBound_(0.0),
    upperBound_(1.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::variableHeightFlowRateFvPatchScalarField
::variableHeightFlowRateFvPatchScalarField
(
    const variableHeightFlowRateFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    lowerBound_(ptf.lowerBound_),
    upperBound_(ptf.upperBound_)
{}


Foam::variableHeightFlowRateFvPatchScalarField
::variableHeightFlowRateFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    lowerBound_(readScalar(dict.lookup("lowerBound"))),
    upperBound_(readScalar(dict.lookup("upperBound")))
{
    patchType() = dict.lookupOrDefault<word>("patchType", word::null);
    this->refValue() = 0.0;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(this->patchInternalField());
    }

    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::variableHeightFlowRateFvPatchScalarField
    ::variableHeightFlowRateFvPatchScalarField
(
    const variableHeightFlowRateFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    phiName_(ptf.phiName_),
    lowerBound_(ptf.lowerBound_),
    upperBound_(ptf.upperBound_)
{}


Foam::variableHeightFlowRateFvPatchScalarField
    ::variableHeightFlowRateFvPatchScalarField
(
    const variableHeightFlowRateFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    lowerBound_(ptf.lowerBound_),
    upperBound_(ptf.upperBound_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::variableHeightFlowRateFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchFieldInDb<surfaceScalarField, scalar>(db(), phiName_);

    scalarField alphap(this->patchInternalField());


    forAll(phip, i)
    {
        if (phip[i] < -SMALL)
        {
            if (alphap[i] < lowerBound_)
            {
                this->refValue()[i] = 0.0;
            }
            else if (alphap[i] > upperBound_)
            {
                this->refValue()[i] = 1.0;
            }
            else
            {
                this->refValue()[i] = alphap[i];
            }

            this->valueFraction()[i] = 1.0;
        }
        else
        {
            this->refValue()[i] = 0.0;
            this->valueFraction()[i] = 0.0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::variableHeightFlowRateFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    if (phiName_ != "phi")
    {
        os.writeEntry("phi", phiName_);
    }
    os.writeEntry("lowerBound", lowerBound_);
    os.writeEntry("upperBound", upperBound_);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        variableHeightFlowRateFvPatchScalarField
    );
}

// ************************************************************************* //
