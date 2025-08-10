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
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/multiphaseGIBalphaInlet/multiphaseGIBalphaInletFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "global/unitConversion/unitConversion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseGIBalphaInletFvPatchScalarField::
multiphaseGIBalphaInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    blendingFactor_(this->size(), 0),
    alphaInlet_(1.0),
    injectDir_(vector::zero),
    angle_(1.0)
{}


Foam::multiphaseGIBalphaInletFvPatchScalarField::
multiphaseGIBalphaInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const scalar& value
)
:
    fvPatchField<scalar>(p, iF, value),
    blendingFactor_(this->size(), 0),
    alphaInlet_(1.0),
    injectDir_(vector::zero),
    angle_(1.0)
{}


Foam::multiphaseGIBalphaInletFvPatchScalarField::
multiphaseGIBalphaInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fvPatchField<scalar>(p, iF, dict, valueRequired),
    blendingFactor_(this->size(), 0),
    alphaInlet_(dict.lookupOrDefault<scalar>("alphaInlet", 1.0)),
    injectDir_(dict.lookup("injectDir")),
    angle_(dict.lookupOrDefault<scalar>("angle", 1.0))
{
    computeBlendingFactor();
}


Foam::multiphaseGIBalphaInletFvPatchScalarField::
multiphaseGIBalphaInletFvPatchScalarField
(
    const multiphaseGIBalphaInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<scalar>(p, iF),
    blendingFactor_(this->size(), 0),
    alphaInlet_(ptf.alphaInlet_),
    injectDir_(ptf.injectDir_),
    angle_(ptf.angle_)
{
    if (ptf.blendingFactor_.size())
    {
        mapper(blendingFactor_, ptf.blendingFactor_);
    }
}


Foam::multiphaseGIBalphaInletFvPatchScalarField::
multiphaseGIBalphaInletFvPatchScalarField
(
    const multiphaseGIBalphaInletFvPatchScalarField& pivpvf
)
:
    fvPatchField<scalar>(pivpvf),
    blendingFactor_(pivpvf.blendingFactor_),
    alphaInlet_(pivpvf.alphaInlet_),
    injectDir_(pivpvf.injectDir_),
    angle_(pivpvf.angle_)
{}


Foam::multiphaseGIBalphaInletFvPatchScalarField::
multiphaseGIBalphaInletFvPatchScalarField
(
    const multiphaseGIBalphaInletFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(pivpvf, iF),
    blendingFactor_(pivpvf.blendingFactor_),
    alphaInlet_(pivpvf.alphaInlet_),
    injectDir_(pivpvf.injectDir_),
    angle_(pivpvf.angle_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseGIBalphaInletFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<scalar>::autoMap(m);

    if (blendingFactor_.size())
    {
        m(blendingFactor_, blendingFactor_);
    }
}


void Foam::multiphaseGIBalphaInletFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fvPatchField<scalar>::rmap(ptf, addr);

    if (blendingFactor_.size())
    {
        const multiphaseGIBalphaInletFvPatchScalarField& tiptf =
            refCast<const multiphaseGIBalphaInletFvPatchScalarField>(ptf);

        blendingFactor_.rmap(tiptf.blendingFactor_, addr);
    }
}


void Foam::multiphaseGIBalphaInletFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fvPatchField<scalar>::autoMapGIB(mapper);
    mapper.map(blendingFactor_, scalar(0));
}


void Foam::multiphaseGIBalphaInletFvPatchScalarField::computeBlendingFactor()
{
    tmp<vectorField> nHat = this->patch().nf();

    blendingFactor_ = pos0(((-injectDir_/mag(injectDir_)) & nHat) - Foam::cos(degToRad(angle_)));
}


void Foam::multiphaseGIBalphaInletFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    computeBlendingFactor();

    scalarField::operator=
    (
        blendingFactor_*alphaInlet_ //+ (1.-blendingFactor_)*(this->patchInternalField())
    );

    fvPatchField<scalar>::evaluate();
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseGIBalphaInletFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    // zero for fixedValue, one for zeroGradient
    scalarField coeffs((1.-blendingFactor_)*0.0);

    return tmp<Field<scalar>>
    (
        new Field<scalar>(coeffs)
    );
}

Foam::tmp<Foam::scalarField>
Foam::multiphaseGIBalphaInletFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    // this for fixedValue, zero for zeroGradient
    return (*this); //*blendingFactor_;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseGIBalphaInletFvPatchScalarField::gradientInternalCoeffs() const
{
    // one*deltaCoeffs for fixedValue, zero for zeroGradient
    return -pTraits<scalar>::one*this->patch().deltaCoeffs(); //*blendingFactor_;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseGIBalphaInletFvPatchScalarField::gradientBoundaryCoeffs() const
{
    // this*deltaCoeffs for fixedValue, zero for zeroGradient
    return this->patch().deltaCoeffs()*(*this); //*blendingFactor_;
}


void Foam::multiphaseGIBalphaInletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("alphaInlet", alphaInlet_);
    os.writeEntry("injectDir", injectDir_);
    os.writeEntry("angle", angle_);
    blendingFactor_.writeEntry("blendingFactor", os);
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::multiphaseGIBalphaInletFvPatchScalarField::operator=
(
    const fvPatchField<scalar>& ptf
)
{
    this->evaluate();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        multiphaseGIBalphaInletFvPatchScalarField
    );
}

// ************************************************************************* //
