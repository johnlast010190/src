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

#include "derivedFvPatchFields/totalPressurePressure/totalPressurePressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "addStaticHead/staticHead.H"


// * * * * * * * * * * * * * Private Functions   * * * * * * * * * * * * * * //

void Foam::totalPressurePressureFvPatchScalarField::computeOtherSources()
{
    if (staticHead_.active())
    {
        pStaticHead_ = staticHead_.computeStaticHead();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalPressurePressureFvPatchScalarField::totalPressurePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    basePressureFvPatchScalarField(p, iF),
    pStaticHead_(),
    staticHead_(*this),
    freeStreamVelocity_(nullptr)
{}


Foam::totalPressurePressureFvPatchScalarField::totalPressurePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    basePressureFvPatchScalarField(p, iF, dict),
    pStaticHead_(),
    staticHead_(*this, p, dict),
    freeStreamVelocity_
    (
        dict.found("freeStreamVelocity")
        ? new point(dict.lookup("freeStreamVelocity"))
        : nullptr
    )
{
    if (staticHead_.active())
    {
        if (dict.found("pStaticHead"))
        {
            pStaticHead_ = scalarField("pStaticHead", dict, p.size());
        }
        else
        {
            pStaticHead_ = scalarField(p.size(), 0);
        }
    }
}


Foam::totalPressurePressureFvPatchScalarField::totalPressurePressureFvPatchScalarField
(
    const totalPressurePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    basePressureFvPatchScalarField(ptf, p, iF, mapper),
    pStaticHead_(mapper(ptf.pStaticHead_)),
    staticHead_(*this, ptf.staticHead_, mapper),

    freeStreamVelocity_
    (
         ptf.freeStreamVelocity_.valid()
       ? new vector(ptf.freeStreamVelocity_())
       : nullptr
    )
{}


Foam::totalPressurePressureFvPatchScalarField::totalPressurePressureFvPatchScalarField
(
    const totalPressurePressureFvPatchScalarField& tppsf
)
:
    basePressureFvPatchScalarField(tppsf),
    pStaticHead_(tppsf.pStaticHead_),
    staticHead_(*this, tppsf.staticHead_),
    freeStreamVelocity_
    (
         tppsf.freeStreamVelocity_.valid()
       ? new vector(tppsf.freeStreamVelocity_())
       : nullptr
    )
{}


Foam::totalPressurePressureFvPatchScalarField::totalPressurePressureFvPatchScalarField
(
    const totalPressurePressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    basePressureFvPatchScalarField(tppsf, iF),
    pStaticHead_(tppsf.pStaticHead_),
    staticHead_(*this, tppsf.staticHead_),
    freeStreamVelocity_
    (
         tppsf.freeStreamVelocity_.valid()
       ? new vector(tppsf.freeStreamVelocity_())
       : nullptr
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Mapping functions

//- Map (and resize as needed) from self given a mapping object
void Foam::totalPressurePressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    basePressureFvPatchScalarField::autoMap(m);

    if (staticHead_.active())
    {
        m(pStaticHead_, pStaticHead_);
        staticHead_.autoMap(m);
    }
}


//- Reverse map the given fvPatchField onto this fvPatchField
void Foam::totalPressurePressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    basePressureFvPatchScalarField::rmap(ptf, addr);

    const totalPressurePressureFvPatchScalarField& tiptf =
        refCast<const totalPressurePressureFvPatchScalarField>(ptf);

    if (staticHead_.active())
    {
        pStaticHead_.rmap(tiptf.pStaticHead_, addr);
        staticHead_.rmap(tiptf.staticHead_, addr);
    }
}


//- Map (and resize as needed) from self given a mapping object
void Foam::totalPressurePressureFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    basePressureFvPatchScalarField::autoMapGIB(mapper);
    if (staticHead_.active())
    {
        mapper.map(pStaticHead_, scalar(0));
        staticHead_.autoMapGIB(mapper);
    }
}



Foam::tmp<Foam::scalarField>
Foam::totalPressurePressureFvPatchScalarField::C0Field() const
{
    tmp<scalarField> tc0(new scalarField(this->p0()));
    scalarField& c0 = tc0.ref();

    if (staticHead_.active())
    {
        c0 += pStaticHead_;
    }

    if (freeStreamVelocity_.valid())
    {
        scalarField c0fs(c0.size(), Zero);

        const vectorField UpIn
        (
            this->db().lookupObject<volVectorField>(UName()).
                boundaryField()[this->patch().index()]
        );
        scalar magSqUfreeSteam = magSqr(freeStreamVelocity_());

        // Urel^2 = (U_i - Ufree_i)^2 = Ui^2-2*U_i&Ufree_i + Ufree_i^2;
        // Ui^2 is added implicitly in c2 and c3
        c0fs = neg(phiP())*((freeStreamVelocity_()&UpIn) - 0.5*magSqUfreeSteam);
        if (internalField().dimensions() == dimPressure)
        {
            c0fs *= this->lookupPatchField<volScalarField, scalar>(rhoName());
        }

        c0 += c0fs;
    }
    return tc0;
}


Foam::tmp<Foam::scalarField>
Foam::totalPressurePressureFvPatchScalarField::C1Field() const
{
    return tmp<scalarField>
    (
        new scalarField(this->size(), Zero)
    );
}


Foam::tmp<Foam::scalarField>
Foam::totalPressurePressureFvPatchScalarField::C2Field() const
{
    tmp<scalarField> tc2
    (
        new scalarField(this->size(), -0.5)
    );
    scalarField& c2 = tc2.ref();
    c2 *= neg(phiP());

    if (internalField().dimensions() == dimPressure)
    {
        c2 *= this->lookupPatchField<volScalarField, scalar>(rhoName());
    }

    return tc2;
}


Foam::tmp<Foam::scalarField>
Foam::totalPressurePressureFvPatchScalarField::C3Field() const
{
    tmp<scalarField> tc3
    (
        new scalarField(this->size(), -0.5)
    );
    scalarField& c3 = tc3.ref();
    c3 *= neg(phiP());

    if (internalField().dimensions() == dimPressure)
    {
        c3 *= this->lookupPatchField<volScalarField, scalar>(rhoName());
    }

    return tc3;
}


void Foam::totalPressurePressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    p0().writeEntry("p0", os);
    writeEntryIfDifferent<word>(os, "U", "U", UName());
    if (!staticHead_.active())
    {
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName());
    }
    else
    {
        pStaticHead_.writeEntry("pStaticHead", os);
        staticHead_.write(os);
    }
    if (freeStreamVelocity_.valid())
    {
        os.writeEntry("freeStreamVelocity", freeStreamVelocity_());
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        totalPressurePressureFvPatchScalarField
    );
}

// ************************************************************************* //
