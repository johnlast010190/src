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

#include "fields/fvPatchFields/derived/surfaceNormalFixedValue/surfaceNormalFixedValueFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_(p.size()),
    profile_(),
    outlet_(true)
{}


Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    refValue_("refValue", dict, p.size()),
    profile_(),
    outlet_(dict.lookupOrDefault<Switch>("implicitOutlet", true))
{
    if (dict.found("profile"))
    {
        profile_ = Function1<scalar>::New("profile", dict);
        const scalar t = this->db().time().timeOutputValue();
        refValue_ = profile_->value(t);
    }

    forceAssign(refValue_*patch().nf());
}


Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const surfaceNormalFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_(ptf.refValue_.size(), pTraits<scalar>::zero),
    profile_(),
    outlet_(ptf.outlet_)
{
    mapper(refValue_, ptf.refValue_);

    if (ptf.profile_.valid())
    {
        profile_.reset(ptf.profile_->clone().ptr());
    }

    // Note: refValue_ will have default value of 0 for unmapped faces. This
    // can temporarily happen during e.g. redistributePar. We should not
    // access ptf.patch() instead since redistributePar has destroyed this
    // at the time of mapping.
    forceAssign(refValue_*patch().nf());
}


Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const surfaceNormalFixedValueFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    refValue_(pivpvf.refValue_),
    profile_(),
    outlet_(pivpvf.outlet_)
{
    if (pivpvf.profile_.valid())
    {
        profile_.reset(pivpvf.profile_->clone().ptr());
    }
}


Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const surfaceNormalFixedValueFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    refValue_(pivpvf.refValue_),
    profile_(),
    outlet_(pivpvf.outlet_)
{
    if (pivpvf.profile_.valid())
    {
        profile_.reset(pivpvf.profile_->clone().ptr());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceNormalFixedValueFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    m(refValue_, refValue_);
}


void Foam::surfaceNormalFixedValueFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
    const surfaceNormalFixedValueFvPatchVectorField& tiptf =
        refCast<const surfaceNormalFixedValueFvPatchVectorField>(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}


void Foam::surfaceNormalFixedValueFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::autoMapGIB(mapper);
    mapper.map(refValue_, scalar(0));
}


void Foam::surfaceNormalFixedValueFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (profile_.valid())
    {
        const scalar t = this->db().time().timeOutputValue();
        refValue_ = profile_->value(t);
    }

    vectorField Up(refValue_*patch().nf());
    frameFieldUpdate(Up);

    fvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::vectorField>
Foam::surfaceNormalFixedValueFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    if (outlet_)
    {
        tmp<vectorField> vic
        (
            new vectorField(this->size(), pTraits<vector>::one)
        );

        vic.ref() *= pos0(refValue_);

        return vic;
    }
    else
    {
        return fixedValueFvPatchField<vector>::valueInternalCoeffs(w);
    }

}

Foam::tmp<Foam::vectorField>
Foam::surfaceNormalFixedValueFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    if (outlet_)
    {
        tmp<vectorField> vbc(new vectorField(*this));
        vbc.ref() *= neg(refValue_);

        return vbc;
    }
    else
    {
        return fixedValueFvPatchField<vector>::valueBoundaryCoeffs(w);
    }
}


Foam::tmp<Foam::vectorField>
Foam::surfaceNormalFixedValueFvPatchVectorField::valueDivInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return fixedValueFvPatchField<vector>::valueInternalCoeffs(w);
}


Foam::tmp<Foam::vectorField>
Foam::surfaceNormalFixedValueFvPatchVectorField::valueDivBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return fixedValueFvPatchField<vector>::valueBoundaryCoeffs(w);
}


void Foam::surfaceNormalFixedValueFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchField<vector>::write(os);
    refValue_.writeEntry("refValue", os);
    if (profile_.valid())
    {
        profile_->writeData(os);
    }
    writeEntryIfDifferent<Switch>(os, "implicitOutlet", true, outlet_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        surfaceNormalFixedValueFvPatchVectorField
    );
}

// ************************************************************************* //
