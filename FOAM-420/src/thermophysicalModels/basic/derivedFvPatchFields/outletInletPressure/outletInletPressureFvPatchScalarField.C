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
    (c) 2017 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/outletInletPressure/outletInletPressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "addStaticHead/staticHead.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outletInletPressureFvPatchScalarField::
outletInletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    outletInletFvPatchScalarField(p, iF),
    outletValueField_(this->size(), 0.0),
    staticHead_(*this)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::outletInletPressureFvPatchScalarField::
outletInletPressureFvPatchScalarField
(
    const outletInletPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    outletInletFvPatchScalarField(ptf, p, iF, mapper),
    outletValueField_(mapper(ptf.outletValueField_)),
    staticHead_(*this, ptf.staticHead_, mapper)
{}


Foam::outletInletPressureFvPatchScalarField::
outletInletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    outletInletFvPatchScalarField(p, iF),
    outletValueField_(scalarField("outletValue", dict, p.size())),
    staticHead_(*this, p, dict)
{
    this->patchType() = dict.lookupOrDefault<word>("patchType", word::null);
    this->refValue() = outletValueField_;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            Field<scalar>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(this->refValue());
    }

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::outletInletPressureFvPatchScalarField::
outletInletPressureFvPatchScalarField
(
    const outletInletPressureFvPatchScalarField& tppsf
)
:
    outletInletFvPatchScalarField(tppsf),
    outletValueField_(tppsf.outletValueField_),
    staticHead_(*this, tppsf.staticHead_)
{}


Foam::outletInletPressureFvPatchScalarField::
outletInletPressureFvPatchScalarField
(
    const outletInletPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    outletInletFvPatchScalarField(tppsf, iF),
    outletValueField_(tppsf.outletValueField_),
    staticHead_(*this, tppsf.staticHead_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::outletInletPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    outletInletFvPatchScalarField::autoMap(m);
    m(outletValueField_, outletValueField_);

    staticHead_.autoMap(m);
}


void Foam::outletInletPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    outletInletFvPatchScalarField::rmap(ptf, addr);

    const outletInletPressureFvPatchScalarField& tiptf =
        refCast<const outletInletPressureFvPatchScalarField>(ptf);

    outletValueField_.rmap(tiptf.outletValueField_, addr);

    staticHead_.rmap(tiptf.staticHead_, addr);
}


void Foam::outletInletPressureFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    outletInletFvPatchScalarField::autoMapGIB(mapper);

    mapper.map(outletValueField_, scalar(0));

    staticHead_.autoMapGIB(mapper);
}


void Foam::outletInletPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    tmp<scalarField> pp ( new scalarField(outletValueField_));

    staticHead_.addStaticHead(pp.ref());

    this->refValue() = pp.ref();

    outletInletFvPatchScalarField::updateCoeffs();
}


void Foam::outletInletPressureFvPatchScalarField::write(Ostream& os)
const
{
    fvPatchField<scalar>::write(os);
    if (phiName_ != "phi")
    {
        os.writeEntry("phi", phiName_);
    }
    outletValueField_.writeEntry("outletValue", os);
    this->writeEntry("value", os);

    //-- for the static head
    staticHead_.write(os);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        outletInletPressureFvPatchScalarField
    );
}

// ************************************************************************* //
