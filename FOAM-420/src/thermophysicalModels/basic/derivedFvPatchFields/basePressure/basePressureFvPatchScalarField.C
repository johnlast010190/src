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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/basePressure/basePressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "derivedFvPatchFields/pressureVelocity/pressureVelocityFvPatchVectorField.H"
#include "cfdTools/general/solutionControl/foamCoupledControl/foamCoupledControl.H"


// * * * * * * * * * * * * * Protected Functions   * * * * * * * * * * * * * //

bool Foam::basePressureFvPatchScalarField::isCoupledSolver() const
{
    if (this->db().foundObject<foamCoupledControl>(solutionControl::typeName))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basePressureFvPatchScalarField::basePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    p0_(p.size(), 0.0),
    UName_("U"),
    rhoName_("rho")
{}


Foam::basePressureFvPatchScalarField::basePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    p0_("p0", dict, p.size()),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(p0_);
    }
}


Foam::basePressureFvPatchScalarField::basePressureFvPatchScalarField
(
    const basePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    p0_(mapper(ptf.p0_)),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


Foam::basePressureFvPatchScalarField::basePressureFvPatchScalarField
(
    const basePressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    p0_(tppsf.p0_),
    UName_(tppsf.UName_),
    rhoName_(tppsf.rhoName_)
{}


Foam::basePressureFvPatchScalarField::basePressureFvPatchScalarField
(
    const basePressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    p0_(tppsf.p0_),
    UName_(tppsf.UName_),
    rhoName_(tppsf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalarField& Foam::basePressureFvPatchScalarField::phiP() const
{
    const fvPatchField<vector>& Up =
        patch().lookupPatchFieldInDb<volVectorField, vector>
        (
            db(),
            UName_
        );
    const pressureVelocityFvPatchVectorField& cUp
        = refCast<const pressureVelocityFvPatchVectorField>(Up);

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
        (
            db(),
            cUp.phiName()
        );
    return phip;
}


void Foam::basePressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(p0_, p0_);
}


void Foam::basePressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const basePressureFvPatchScalarField& tiptf =
        refCast<const basePressureFvPatchScalarField>(ptf);
    p0_.rmap(tiptf.p0_, addr);
}


void Foam::basePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField Up
    (
        this->db().lookupObject<volVectorField>(UName_).
            boundaryField()[this->patch().index()]
    );
    if (coorFramePtr())
    {
        makeRelative(Up);
    }

    const vectorField nf(this->patch().nf()());
    const scalarField Un(Up&nf);
    const vectorField Ut(Up - Un*nf);

    // Compute quantities (static head etc) based on the derived class
    computeOtherSources();

    scalarField dp(this->size(), 0.0);
    dp +=
        C0Field()
      + C1Field()*Un
      + C2Field()*magSqr(Un)
      + C3Field()*magSqr(Ut);

    if (coorFramePtr())
    {
        // Relaxation coefficient has to be investigated
        const scalar relax = 0.1;
        forceAssign(*this + relax*(dp - *this));
    }
    else
    {
        forceAssign(dp);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::basePressureFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>& vbc
) const
{
    if (isCoupledSolver() && !coorFramePtr())
    {
        return tmp<Field<scalar>>
        (
            new Field<scalar>(this->size(), Zero)
        );
    }
    else
    {
        return fixedValueFvPatchScalarField::valueBoundaryCoeffs(vbc);
    }
}


Foam::tmp<Foam::scalarField>
Foam::basePressureFvPatchScalarField::gradientBoundaryCoeffs() const
{
    if (isCoupledSolver() && !coorFramePtr())
    {
        return tmp<scalarField>
        (
            new scalarField(this->size(), 0.0)
        );
    }
    else
    {
        return fixedValueFvPatchScalarField::gradientBoundaryCoeffs();
    }
}


void Foam::basePressureFvPatchScalarField::write(Ostream& os) const
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        basePressureFvPatchScalarField
    );
}

// ************************************************************************* //
